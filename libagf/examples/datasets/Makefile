SHELL=/bin/bash

#number of classes in each dataset:
NCLS_heart=2
NCLS_shuttle=7
NCLS_covtype=7
NCLS_segment=7
NCLS_sat=6
NCLS_vehicle=4
NCLS_poker=10
NCLS_mnist=10
NCLS_letter=26
NCLS_codrna=2
NCLS_protein=3
NCLS_pendigits=10
NCLS_dna=3
NCLS_usps=10
NCLS_seismic=3
NCLS_ijcnn1=2
NCLS_leu=2
NCLS_mushrooms=2
NCLS_madelon=2
NCLS_phishing=2
NCLS_splice=2
#this one actually sets the number of classes:
NCLS_humidity=8

NCLS=$(NCLS_$(VER))

#directory containing unprocessed data:
DATADIR=../statlog/data/$(VER)

#working directory for processing data:
WORKDIR=.

#file names:
BASE=$(WORKDIR)/$(VER)
#TRAIN=$(BASE).trn
#TEST=$(BASE).tst

#process training data:
$(BASE).vec $(BASE).cls: $(BASE).0.cls $(BASE).0.vec
	agf_preprocess $(PRE) -U $(BASE).0 $(BASE)

#convert to LIBSVM format:
$(BASE).svm: $(BASE).cls $(BASE).vec
	agf2ascii -M $(BASE) $(BASE).svm

#file conversion:
#steps:
#1. unify test and training, if applicable
#2. convert to libAGF format and copy to working directory
#3. pre-processing: 
#     - normalization (optional)
#     - SVD (optional)
#     - make class labels go from 0 to Nc-1
#4. convert to LIBSVM format (if needed) above


#list of datasets from LIBSVM website with test and training:
SVMTTDATA=madelon \
	ijcnn1 \
	seismic \
	usps \
	dna \
	pendigits \
	codrna \
	poker 
	#splice

#list of datasets that need to be separated into test and training:
ALLINONE=	heart \
		segment \
		mushrooms \
		phishing \
		covtype \
		letter \
		humidity

#large datasets that we want to subsample (because we haven't got all year...):
$(WORKDIR)/mnist_s.0.cls:
	make VER=mnist
	agf_preprocess -zf 0.1 $(WORKDIR)/mnist.0 dum $(WORKDIR)/mnist_s.0

#shitty way of doing it...
$(WORKDIR)/covtype_s.0.cls:
	make VER=covtype
	agf_preprocess -zf 0.01 $(WORKDIR)/covtype.0 dum $(WORKDIR)/covtype_s.0

#make the working directory:
#$(WORKDIR)/$(VER).timestamp:
#	if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
#	date > $(WORKDIR)/$(VER).timestamp

#UCI datasets in LVQ format (sans header) with only training data:
$(WORKDIR)/heart.0.cls: $(WORKDIR)/%.cls: $(DATADIR)/%.dat.txt 
	#$(WORKDIR)/$(VER).timestamp
	lvq2agf -H $< $(WORKDIR)/$(VER).0

$(WORKDIR)/segment.0.cls: $(DATADIR)/segment.dat.txt 
	#$(WORKDIR)/$(VER).timestamp
	lvq2agf -H $< $(WORKDIR)/$(VER).0
	#3rd feature must be removed since it only takes on one value:
	agf_preprocess -F $(BASE).0 $(BASE).0 0 1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18

$(WORKDIR)/letter.0.cls: $(DATADIR)/letter-recognition.data.txt 
	#$(WORKDIR)/$(VER).timestamp
	sed 's/A/0/;s/B/1/;s/C/2/;s/D/3/;s/E/4/;s/F/5/;s/G/6/;s/H/7/;s/I/8/;s/J/9/;s/K/10/;s/L/11/;s/M/12/;s/N/13/;s/O/14/;s/P/15/;s/Q/16/;s/Q/16/;s/R/17/;s/S/18/;s/T/19/;s/U/20/;s/V/21/;s/W/22/;s/X/23/;s/Y/24/;s/Z/25/;' < $< | awk -F ',' '{print $$2, $$3, $$4, $$5, $$6, $$7, $$8, $$9, $$10, $$11, $$12, $$13, $$14, $$15, $$16, $$1}' > $(WORKDIR)/letter.txt
	lvq2agf -H $(WORKDIR)/letter.txt $(WORKDIR)/letter.0

#UCI datasets in LVQ format (sans header) with both training and test data:
$(WORKDIR)/sat.0.cls $(WORKDIR)/shuttle.0.cls: $(WORKDIR)/%.0.cls: $(DATADIR)/%.trn.txt $(DATADIR)/%.tst.txt
	#$(WORKDIR)/$(VER).timestamp
	cat $(DATADIR)/$(VER).trn.txt $(DATADIR)/$(VER).tst.txt > $(BASE).0.txt
	lvq2agf -H $(BASE).0.txt $(BASE).0

#$(WORKDIR)/sat.tst.raw.cls $(WORKDIR)/shuttle.tst.raw.cls: $(WORKDIR)/%.tst.raw.cls: $(DATADIR)/%.tst.txt 
	#$(WORKDIR)/$(VER).timestamp
#	lvq2agf -H $(DATADIR)/$(VER).tst.txt $(TEST).raw

#LIBSVM datasets with only training data:
$(WORKDIR)/mushrooms.0.cls $(WORKDIR)/phishing.0.cls $(WORKDIR)/covtype.0.cls: $(WORKDIR)/%.0.cls: $(DATADIR)/% 
	#$(WORKDIR)/$(VER).timestamp
	svm2agf -E 0 $(DATADIR)/$(VER) $(BASE).0

#LIBSVM datasets with both test and training data:
$(addsuffix .0.cls, $(addprefix $(WORKDIR)/, $(SVMTTDATA))): $(WORKDIR)/%.0.cls: $(DATADIR)/% $(DATADIR)/%.t
	#$(WORKDIR)/$(VER).timestamp
	cat $(DATADIR)/$(VER) $(DATADIR)/$(VER).t > $(BASE).0.svm
	svm2agf -E 0 $(BASE).0.svm $(BASE).0

#$(addsuffix .tst.raw.cls, $(addprefix $(WORKDIR)/, $(SVMTTDATA))): $(WORKDIR)/%.tst.raw.cls: $(DATADIR)/%.t 
	#$(WORKDIR)/$(VER).timestamp
#	svm2agf -E 0 $(DATADIR)/$(VER).t $(TEST).raw

#rearrange splice features:
SED_PROGRAM='s/1\./3\./;s/2\./2\./;s/3\./0\./;s/4\./1\./;'
$(WORKDIR)/splice.0.cls: $(DATADIR)/splice $(DATADIR)/splice.t
	#$(WORKDIR)/$(VER).timestamp
	cat $(DATADIR)/slice $(DATADIR)/splice.t > $(BASE).0.svm
	sed $(SED_PROGRAM) < $(BASE).0.svm > $(BASE).0.svm
	svm2agf -E 0 $(WORKDIR)/splice.txt $(BASE).0

#$(WORKDIR)/splice.0.cls: $(DATADIR)/splice.t 
	#$(WORKDIR)/$(VER).timestamp
#	sed $(SED_PROGRAM) < $< > $(WORKDIR)/splice.t.txt
#	svm2agf -E 0 $(WORKDIR)/splice.t.txt $(WORKDIR)/splice.tst.raw

$(WORKDIR)/mnist.0.cls: $(DATADIR)/mnist $(DATADIR)/mnist.t 
	#$(WORKDIR)/$(VER).timestamp
	cat $(DATADIR)/$(VER) $(DATADIR)/$(VER).t > $(BASE).0.svm
	svm2agf -E 0 $(BASE).0.svm $(BASE).0
	agf_preprocess -F $(BASE).0 $(BASE).0 $(MNIST_INDICES)

#my own dataset:
$(WORKDIR)/humidity.0.cls: ../humidity_data/train_SQ36s.dat ../humidity_data/est_bt_cld_land_z0.00.vec
	cp ../humidity_data/est_bt_cld_land_z0.00.vec $(WORKDIR)/humidity.0.vec
	float_to_class -q $(NCLS_humidity) -g ../humidity_data/train_SQ36s.dat $(WORKDIR)/humidity.0.cls 0.00007 0.001

$(WORKDIR)/vehicle.0.cls: $(DATADIR)/xa?.dat.txt
	cat $(DATADIR)/xa?.dat.txt | sed 's/opel/0/;s/saab/1/;s/bus/2/;s/van/3/' > $(BASE).0.txt
	lvq2agf -H $(BASE).0.txt $(BASE).0

#we like the mnist dataset because it takes 3 days to train in LIBSVM:
MNIST_INDICES = 33 \
34 \
35 \
36 \
37 \
38 \
39 \
40 \
41 \
42 \
43 \
44 \
45 \
46 \
47 \
48 \
61 \
62 \
63 \
64 \
65 \
66 \
67 \
68 \
69 \
70 \
71 \
72 \
73 \
74 \
75 \
76 \
77 \
78 \
79 \
80 \
88 \
89 \
90 \
91 \
92 \
93 \
94 \
95 \
96 \
97 \
98 \
99 \
100 \
101 \
102 \
103 \
104 \
105 \
106 \
107 \
108 \
109 \
110 \
115 \
116 \
117 \
118 \
119 \
120 \
121 \
122 \
123 \
124 \
125 \
126 \
127 \
128 \
129 \
130 \
131 \
132 \
133 \
134 \
135 \
136 \
137 \
138 \
142 \
143 \
144 \
145 \
146 \
147 \
148 \
149 \
150 \
151 \
152 \
153 \
154 \
155 \
156 \
157 \
158 \
159 \
160 \
161 \
162 \
163 \
164 \
165 \
166 \
170 \
171 \
172 \
173 \
174 \
175 \
176 \
177 \
178 \
179 \
180 \
181 \
182 \
183 \
184 \
185 \
186 \
187 \
188 \
189 \
190 \
191 \
192 \
193 \
194 \
195 \
197 \
198 \
199 \
200 \
201 \
202 \
203 \
204 \
205 \
206 \
207 \
208 \
209 \
210 \
211 \
212 \
213 \
214 \
215 \
216 \
217 \
218 \
219 \
220 \
221 \
222 \
223 \
225 \
226 \
227 \
228 \
229 \
230 \
231 \
232 \
233 \
234 \
235 \
236 \
237 \
238 \
239 \
240 \
241 \
242 \
243 \
244 \
245 \
246 \
247 \
248 \
249 \
250 \
251 \
252 \
253 \
254 \
255 \
256 \
257 \
258 \
259 \
260 \
261 \
262 \
263 \
264 \
265 \
266 \
267 \
268 \
269 \
270 \
271 \
272 \
273 \
274 \
275 \
276 \
277 \
278 \
279 \
280 \
281 \
282 \
283 \
284 \
285 \
286 \
287 \
288 \
289 \
290 \
291 \
292 \
293 \
294 \
295 \
296 \
297 \
298 \
299 \
300 \
301 \
302 \
303 \
304 \
305 \
306 \
307 \
308 \
309 \
310 \
311 \
312 \
313 \
314 \
315 \
316 \
317 \
318 \
319 \
320 \
321 \
322 \
323 \
324 \
325 \
326 \
327 \
328 \
329 \
330 \
331 \
332 \
333 \
334 \
335 \
337 \
338 \
339 \
340 \
341 \
342 \
343 \
344 \
345 \
346 \
347 \
348 \
349 \
350 \
351 \
352 \
353 \
354 \
355 \
356 \
357 \
358 \
359 \
360 \
361 \
362 \
365 \
366 \
367 \
368 \
369 \
370 \
371 \
372 \
373 \
374 \
375 \
376 \
377 \
378 \
379 \
380 \
381 \
382 \
383 \
384 \
385 \
386 \
387 \
388 \
389 \
390 \
394 \
395 \
396 \
397 \
398 \
399 \
400 \
401 \
402 \
403 \
404 \
405 \
406 \
407 \
408 \
409 \
410 \
411 \
412 \
413 \
414 \
415 \
416 \
417 \
418 \
422 \
423 \
424 \
425 \
426 \
427 \
428 \
429 \
430 \
431 \
432 \
433 \
434 \
435 \
436 \
437 \
438 \
439 \
440 \
441 \
442 \
443 \
444 \
445 \
446 \
447 \
449 \
450 \
451 \
452 \
453 \
454 \
455 \
456 \
457 \
458 \
459 \
460 \
461 \
462 \
463 \
464 \
465 \
466 \
467 \
468 \
469 \
470 \
471 \
472 \
473 \
474 \
475 \
478 \
479 \
480 \
481 \
482 \
483 \
484 \
485 \
486 \
487 \
488 \
489 \
490 \
491 \
492 \
493 \
494 \
495 \
496 \
497 \
498 \
499 \
500 \
501 \
502 \
503 \
505 \
506 \
507 \
508 \
509 \
510 \
511 \
512 \
513 \
514 \
515 \
516 \
517 \
518 \
519 \
520 \
521 \
522 \
523 \
524 \
525 \
526 \
527 \
528 \
529 \
530 \
531 \
533 \
534 \
535 \
536 \
537 \
538 \
539 \
540 \
541 \
542 \
543 \
544 \
545 \
546 \
547 \
548 \
549 \
550 \
551 \
552 \
553 \
554 \
555 \
556 \
557 \
558 \
559 \
561 \
562 \
563 \
564 \
565 \
566 \
567 \
568 \
569 \
570 \
571 \
572 \
573 \
574 \
575 \
576 \
577 \
578 \
579 \
580 \
581 \
582 \
583 \
584 \
585 \
586 \
590 \
591 \
592 \
593 \
594 \
595 \
596 \
597 \
598 \
599 \
600 \
601 \
602 \
603 \
604 \
606 \
607 \
608 \
609 \
610 \
611 \
612 \
613 \
614 \
618 \
619 \
620 \
621 \
622 \
623 \
624 \
625 \
626 \
627 \
628 \
629 \
630 \
631 \
632 \
633 \
634 \
635 \
636 \
637 \
638 \
639 \
640 \
641 \
642 \
646 \
647 \
648 \
649 \
650 \
651 \
652 \
653 \
654 \
655 \
656 \
657 \
658 \
659 \
660 \
661 \
662 \
663 \
664 \
665 \
666 \
667 \
668 \
669 \
670 \
675 \
676 \
677 \
678 \
679 \
680 \
681 \
682 \
683 \
684 \
685 \
686 \
687 \
688 \
689 \
690 \
691 \
692 \
693 \
694 \
695 \
696 \
697 \
703 \
704 \
705 \
706 \
707 \
708 \
709 \
710 \
711 \
712 \
713 \
714 \
715 \
716 \
717 \
718 \
719 \
720 \
721 \
722 \
723 \
724 \
725 \
732 \
733 \
734 \
735 \
736 \
737 \
738 \
739 \
740 \
741 \
742 \
743 \
744 \
745 \
746 \
747 \
748 \
749 \
750 \
751 \
752 \
762 \
763 \
764 \
765 \
766 \
767 \
768 \
769 \
770 \
771 \
772 \
773 \
774 \
775 \
776 \
777 

