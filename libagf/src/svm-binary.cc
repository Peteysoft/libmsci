#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*
 * #include <ctype.h>
#include <float.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <locale.h>
#include <errno.h>
*/

#include "svm.h"

#ifndef min
template <class T> static inline T min(T x,T y) { return (x<y)?x:y; }
#endif
#ifndef max
template <class T> static inline T max(T x,T y) { return (x>y)?x:y; }
#endif

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

static inline double powi(double base, int times)
{
	double tmp = base, ret = 1.0;

	for(int t=times; t>0; t/=2)
	{
		if(t%2==1) ret*=tmp;
		tmp = tmp * tmp;
	}
	return ret;
}

//returns the dot product:
double svm_dot(const svm_node *px, const svm_node *py)
{
	double sum = 0;
	while(px->index != -1 && py->index != -1)
	{
		if(px->index == py->index)
		{
			sum += px->value * py->value;
			++px;
			++py;
		}
		else
		{
			if(px->index > py->index)
				++py;
			else
				++px;
		}			
	}
	return sum;
}

double svm_distance2(const svm_node *x, const svm_node *y, svm_node *deriv) {
	double sum = 0;
	while(x->index != -1 && y->index !=-1)
	{
		if(x->index == y->index)
		{
			double d = x->value - y->value;
			sum += d*d;
			deriv->index=x->index;
			deriv->value=2*d;
			++x;
			++y;
			++deriv;
		}
		else
		{
			if(x->index > y->index)
			{	
				sum += y->value * y->value;
				deriv->index=y->index;
				deriv->value=2*y->value;
				++y;
				++deriv;
			}
			else
			{
				sum += x->value * x->value;
				deriv->index=x->index;
				deriv->value=-2*x->value;
				++x;
				++deriv;
			}
		}
	}

	while(x->index != -1)
	{
		sum += x->value * x->value;
		deriv->index=x->index;
		deriv->value=-2*x->value;
		++x;
		++deriv;
	}

	while(y->index != -1)
	{
		sum += y->value * y->value;
		deriv->index=y->index;
		deriv->value=-2*y->value;
		++deriv;
		++y;
	}

	return sum;
}
			
//returns the derivative dot product wrt the first argument:
void svm_dot_deriv(const svm_node *px, const svm_node *py, svm_node *deriv)
{
	double sum = 0;
	while(px->index != -1 && py->index != -1)
	{
		if(px->index == py->index)
		{
			deriv->index=px->index;
			deriv->value=py->value;
			++px;
			++py;
			++deriv;
		}
		else
		{
		//if derivative lacks index, just assume it's zero!
			if(px->index > py->index)
				++py;
			else 
				++px;
		}			
	}
	deriv->index=-1;
}

double k_function_deriv(const svm_node *x, const svm_node *y,
			  const svm_parameter& param, svm_node *deriv)
{
	switch(param.kernel_type)
	{
		case LINEAR:
			svm_dot_deriv(x, y, deriv);
			return svm_dot(x,y);
		case POLY:
			double t1=param.gamma*svm_dot(x,y)+param.coef0;
			double t2=powi(t1, param.degree-1)*param.degree;
			svm_dot_deriv(x, y, deriv);
			while (deriv->index!=-1) {
				deriv->value*=t2;
				++deriv;
			}
			return powi(t1, param.degree);
		case RBF:
		{
			double t1 = exp(-param.gamma*svm_distance2(x, y, deriv));
			while (deriv->index!=-1) {
				deriv->value*=-param.gamma*t1;
				++deriv;
			}
			return t2;
		}
		case SIGMOID:
			double t1=tanh(param.gamma*dot(x,y)+param.coef0);
			double t2=param.gamma*(1-t1*t1);
			svm_dot_deriv(x, y, deriv);
			while (deriv->index!=-1) {
				deriv->value*=t2;
				++deriv;
			}
			return t1;
		case PRECOMPUTED:  //x: test (validation), y: SV
			fprintf(stderr, "Pre-computed kernels not supported!\n");
			exit(1);
		default:
			return 0;  // Unreachable 
	}
}

double svm_predict_deriv(const svm_model *model, const svm_node *x, svm_node *deriv)
{
	int i;
	int nr_class = model->nr_class;
	double result;
	if(model->param.svm_type == ONE_CLASS ||
	   model->param.svm_type == EPSILON_SVR ||
	   model->param.svm_type == NU_SVR ||
	   nr_class!=2)
	{
		fprintf(stderr, "Only valid for binary classifications\n");
		exit(1);
	}
	else
	{
		int l = model->l;
		double sum=0;
		//yes, we will actually initialize the derivative vector to zero for all indices all
		//the way up to the max allowed--that's how brain-dead this system is...:
		for (i=0; i<max_nr_attr; i++) {
			deriv[i].index=i;
			deriv[i].value=0;
		}
		for(i=0;i<l;i++) {
			double kvalue;
			svm_node deriv1;

			kvalue = k_function_deriv(x,model->SV[i], model->param, &deriv1);
			sum+=model->sv_coef[0][i];
			k=0;
			for (j=0; deriv1[j].index!=-1; j++) {
				if (deriv[k].index==deriv1[j].index) {
					deriv[k].value+=model->sv_coef[0][i]*deriv1[j].value;
					k++;
				} else {
					for (;deriv[k].index<deriv1[k].index; k++);
				}
			}
		}
		sum-=model->rho[p];

		return sum;
	}
}

static double sigmoid_predict(double decision_value, double A, double B, double *deriv)
{
	double expfApB = exp(decision_value*A+B);
	*deriv=-A*expfApB/(1+expfApB)/(1+expfApB);
	return 1.0/(1+exp(fApB)) ;
}

double svm_predict_binary(
	const svm_model *model, const svm_node *x, double *prob_estimates)
{
	int nr_class = model->nr_class;
	assert(nr_class==2);

	if ((model->param.svm_type == C_SVC || model->param.svm_type == NU_SVC) &&
	    model->probA!=NULL && model->probB!=NULL)
	{
		int i;

		double *dec_values = Malloc(double, nr_class*(nr_class-1)/2);
		svm_predict_values(model, x, dec_values);

		double min_prob=1e-7;
		double **pairwise_prob=Malloc(double *,nr_class);
		for(i=0;i<nr_class;i++)
			pairwise_prob[i]=Malloc(double,nr_class);
		int k=0;
		for(i=0;i<nr_class;i++)
			for(int j=i+1;j<nr_class;j++)
			{
				pairwise_prob[i][j]=min(max(sigmoid_predict(dec_values[k],model->probA[k],model->probB[k]),min_prob),1-min_prob);
				pairwise_prob[j][i]=1-pairwise_prob[i][j];
				k++;
			}
		multiclass_probability(nr_class,pairwise_prob,prob_estimates);

		int prob_max_idx = 0;
		for(i=1;i<nr_class;i++)
			if(prob_estimates[i] > prob_estimates[prob_max_idx])
				prob_max_idx = i;
		for(i=0;i<nr_class;i++)
			free(pairwise_prob[i]);
		free(dec_values);
		free(pairwise_prob);
		return model->label[prob_max_idx];
	}
	else {
		fprintf(stderr, "Binary classification models only supported\n");
		exit(1);
	}
}

//main part starts here...

int print_null(const char *s,...) {return 0;}

static int (*info)(const char *fmt,...) = &printf;

struct svm_node *x;
int max_nr_attr = 64;

struct svm_model* model;
int predict_probability=1;

static char *line = NULL;
static int max_line_len;

//why isn't this in the main library??
static char* readline(FILE *input)
{
	int len;

	if(fgets(line,max_line_len,input) == NULL)
		return NULL;

	while(strrchr(line,'\n') == NULL)
	{
		max_line_len *= 2;
		line = (char *) realloc(line,max_line_len);
		len = (int) strlen(line);
		if(fgets(line+len,max_line_len-len,input) == NULL)
			break;
	}
	return line;
}

void exit_input_error(int line_num)
{
	fprintf(stderr,"Wrong input format at line %d\n", line_num);
	exit(1);
}

void predict_binary(FILE *input, FILE *output)
{
	int correct = 0;
	int total = 0;
	double error = 0;
	double sump = 0, sumt = 0, sumpp = 0, sumtt = 0, sumpt = 0;

	int svm_type=svm_get_svm_type(model);
	int nr_class=svm_get_nr_class(model);
	double *prob_estimates=NULL;
	int j;

	if (svm_type==NU_SVR || svm_type==EPSILON_SVR || nr_class!=2 || predict_probability!=1)
	{
		info("Binary class with probabilities estimates only supported.\n");
		exit(1);
	}
	else
	{
		int *labels=(int *) malloc(nr_class*sizeof(int));
		svm_get_labels(model,labels);
		prob_estimates = (double *) malloc(nr_class*sizeof(double));
		fprintf(output,"labels");		
		for(j=0;j<nr_class;j++)
			fprintf(output," %d",labels[j]);
		fprintf(output,"\n");
		free(labels);
	}

	max_line_len = 1024;
	line = (char *)malloc(max_line_len*sizeof(char));
	while(readline(input) != NULL)
	{
		int i = 0;
		double target_label, predict_label;
		char *idx, *val, *label, *endptr;
		int inst_max_index = -1; // strtol gives 0 if wrong format, and precomputed kernel has <index> start from 0

		label = strtok(line," \t\n");
		if(label == NULL) // empty line
			exit_input_error(total+1);

		target_label = strtod(label,&endptr);
		if(endptr == label || *endptr != '\0')
			exit_input_error(total+1);

		while(1)
		{
			if(i>=max_nr_attr-1)	// need one more for index = -1
			{
				max_nr_attr *= 2;
				x = (struct svm_node *) realloc(x,max_nr_attr*sizeof(struct svm_node));
			}

			idx = strtok(NULL,":");
			val = strtok(NULL," \t");

			if(val == NULL)
				break;
			errno = 0;
			x[i].index = (int) strtol(idx,&endptr,10);
			if(endptr == idx || errno != 0 || *endptr != '\0' || x[i].index <= inst_max_index)
				exit_input_error(total+1);
			else
				inst_max_index = x[i].index;

			errno = 0;
			x[i].value = strtod(val,&endptr);
			if(endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
				exit_input_error(total+1);

			++i;
		}
		x[i].index = -1;

		predict_label = svm_predict_binary(model,x,prob_estimates);
		fprintf(output,"%g",predict_label);
		for(j=0;j<nr_class;j++)
			fprintf(output," %.16lg",prob_estimates[j]);
		fprintf(output,"\n");

		if(predict_label == target_label)
			++correct;
		error += (predict_label-target_label)*(predict_label-target_label);
		sump += predict_label;
		sumt += target_label;
		sumpp += predict_label*predict_label;
		sumtt += target_label*target_label;
		sumpt += predict_label*target_label;
		++total;
	}
	info("Accuracy = %g%% (%d/%d) (classification)\n",
			(double)correct/total*100,correct,total);
	if(predict_probability)	free(prob_estimates);
}

void exit_with_help()
{
	printf(
	"Usage: svm-binary [options] test_file model_file output_file\n"
	"options:\n"
	"-q : quiet mode (no outputs)\n"
	);
	exit(1);
}

int main(int argc, char **argv)
{
	FILE *input, *output;
	int i;
	predict_probability=1;
	// parse options
	for(i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break;
		++i;
		switch(argv[i-1][1])
		{
			case 'q':
				info = &print_null;
				i--;
				break;
			default:
				fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
		}
	}

	if(i>=argc-2)
		exit_with_help();

	input = fopen(argv[i],"r");
	if(input == NULL)
	{
		fprintf(stderr,"can't open input file %s\n",argv[i]);
		exit(1);
	}

	output = fopen(argv[i+2],"w");
	if(output == NULL)
	{
		fprintf(stderr,"can't open output file %s\n",argv[i+2]);
		exit(1);
	}

	if((model=svm_load_model(argv[i+1]))==0)
	{
		fprintf(stderr,"can't open model file %s\n",argv[i+1]);
		exit(1);
	}

	x = (struct svm_node *) malloc(max_nr_attr*sizeof(struct svm_node));
	if(svm_check_probability_model(model)==0)
	{
		fprintf(stderr,"Model does not support probabiliy estimates\n");
		exit(1);
	}

	predict(input,output);
	svm_free_and_destroy_model(&model);
	free(x);
	free(line);
	fclose(input);
	fclose(output);
	return 0;
}
