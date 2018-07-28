#!/bin/bash

modfile=$1
testfile=$2

test1=svmtest.$RANDOM
test2=svmtest.$RANDOM

#if classify_s has been compiled with double-precision calculations
#the results in both cases should be identical

#solve for conditional probabilities:
svm-predict -b 1 $testfile $modfile $test1.svmout
svmout2agf $test1.svmout $test1
classify_s -ZMA $modfile $testfile $test2.svmout
svmout2agf $test2.svmout $test2
cls_comp_stats $test1.cls $test2

#voting solution, no probabilities returned:
svm-predict -b 0 $testfile $modfile $test1.svmout
svmout2agf $test1.svmout $test1
classify_s -ZMA -Q 1 $modfile $testfile $test2.svmout
svmout2agf $test2.svmout $test2
cls_comp_stats $test1.cls $test2

rm $test1.*
rm $test2.*

