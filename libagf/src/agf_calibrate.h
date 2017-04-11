#ifndef __AGF_CALIBRATE__H__
#define __AGF_CALIBRATE__H__

namespace libagf {
  //uncertainty coefficient for a binary classifier:
  template <typename real>
  real uc_bin(real ntn,			//number of true negatives
		real nfp,		// " false positives
		real nfn,		// " false negatives
		real ntp);		// " true positives

  //accuracy of a binary classifier:
  template <typename real>
  real acc_bin(real ntn,		//number of true negatives
		real nfp,		// " false positives
		real nfn,		// " false negatives
		real ntp);		// " true positives

  //correlation for a binary classifier:
  template <typename real>
  real corr_bin(real ntn,		//number of true negatives
		real nfp,		// " false positives
		real nfn,		// " false negatives
		real ntp);		// " true positives

  //optimizes the skill of a binary classifier using a golden-section search
  template <typename real>
  real optimize_binary_skill(nel_ta *nonecum, 	//cumulative class values
		  real *rsort, 			//sorted probabilities
		  nel_ta n,			//number of results
		  //function that computes skill for binary classifier:
		  real (*binary_skill) (real, real, real, real),
		  real tol, 			//desired tolerance
		  long maxiter);		//maximum number of iterations

  //optimizes the skill of a binary classifier rigourously, by going through the probabilities 
  //one-by-one using each as the threshold value and calculating the skill score
  template <typename real>
  real optimize_binary_skill_rig(nel_ta *nonecum,//cumulative class values
		  real *rsort, 			//sorted probabilities
		  nel_ta n,			//number of results
		  //function that computes skill for binary classifier:
		  real (*binary_skill) (real, real, real, real));

  //optimizes the skill of a binary classifier rigourously, by going through the probabilities 
  //one-by-one using each as the threshold value and calculating the skill score
  template <typename real>
  real *calc_skill_vs_r0(nel_ta *nonecum, 
		  real *rsort,
		  nel_ta n, 
		  real (*binary_skill) (real, real, real, real));

  //sorts the conditional probabilities, cumulates the class labels:
  template <typename real, typename cls_t>
  void sortr_cumulate_ones(cls_t *truth, 	//true class values {0,1}
		  real *r,			//calculated probabilities P(1|x)-P(0|x)
		  nel_ta n,			//number of results
		  nel_ta *nonecum,		//returned cumulated class labels
		  real *rsort);			//returned sorted probabilities

  //integrate the ROC curve:
  template <typename real>
  real integrate_roc(nel_ta *nonecum, 		//cumulated class labels
		  real *rsort, 			//sorted probabilities
		  nel_ta n);			//number of results

  //calculate the ROC curve:
  //(returns results in a two-dimensional array
  //first row is false-alarm-rate, second is hit-rate)
  template <typename real>
  real ** calc_roc(nel_ta *nonecum, 		//cumulated class labels
		  real *rsort,			//sorted probabilities
		  nel_ta n);			//number of results

} //end namespace libagf

#endif

