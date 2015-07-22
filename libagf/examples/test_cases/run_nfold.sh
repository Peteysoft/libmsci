set -e

cp ../humidity_data/train_SQ36s.dat ../humidity_data/est_bt_cld_land_z0.00.dat

agf_preprocess -d 10 -L ../humidity_data/est_bt_cld_land_z0.00 foo
float_to_class -q 8 -g foo-00.dat foo-00.cls 0.00007 0.001
agf_preprocess -d 10 ../humidity_data/est_bt_cld_land_z0.00 bar
agf2ascii foo-00 > check1.txt
agf2ascii bar-00 > check2.txt
diff check1.txt check2.txt
echo "CASE1a" > results.txt

agf_preprocess -f 0.1 -L ../humidity_data/est_bt_cld_land_z0.00 foo bar
float_to_class -q 8 -g bar.dat bar.cls 0.00007 0.001
agf_preprocess -f 0.1 ../humidity_data/est_bt_cld_land_z0.00 foo2 bar2
agf2ascii bar > check1.txt
agf2ascii bar2 > check2.txt
diff check1.txt check2.txt
echo "CASE1b" >> results.txt

multi_borders -s 1 test_case1.txt foo bar foobar.txt > foo_script.sh
diff test_case1.out.sh foo_script.sh >> results.txt
diff test_case1.out.txt foobar.txt >> results.txt
echo "CASE1c" >> results.txt

nfold ../sample_classes/sctrain >> results.txt
echo "CASE2a" >> results.txt
nfold -c 1 -k 200 -k 25 ../sample_classes/sctrain >> results.txt
echo "CASE2b" >> results.txt
nfold -c 2 ../sample_classes/sctrain >> results.txt
echo "CASE2c" >> results.txt
nfold -n -W 50 -k 400 ../sample_classes/sctrain >> results.txt
echo "CASE2d" >> results.txt
nfold -c 2 -k 19 ../sample_classes/sctrain >> results.txt
echo "CASE2e" >> results.txt

agf_preprocess -z -f 0.1 ../humidity_data/est_bt_cld_land_z0.00 dum sq_reduced
nfold -S 5 -c 5 -d 4 ../humidity_data/sq_control4.txt sq_reduced >> results.txt
echo "CASE3a" >> results.txt
nfold -n -S 5 -c 1 -d 4 -k 400 -W 40 sq_reduced >> results.txt
echo "CASE3b" >> results.txt
nfold -c 2 -d 4 -k 3 sq_reduced >> results.txt
echo "CASE3c" >> results.txt
nfold -n -u -c 5 -d 4 ../humidity_data/sq_control4.txt sq_reduced >> results.txt
echo "CASE3d" >> results.txt

cluster_knn -p 2.5 ../sample_classes/sctest.vec foo.cls >> results.txt
cls_comp_stats ../sample_classes/sctest.cls foo >> results.txt
echo "CASE4" >> results.txt

validate_pdf.sh ../sample_classes/sctest.vec >> results.txt
echo "CASE5a" >> results.txt

validate_pdf.sh -k 400 -W 50 ../sample_classes/sctest.vec >> results.txt
echo "CASE5b" >> results.txt

validate_pdf.sh -c 7 ../sample_classes/sctest.vec >> results.txt
echo "CASE5c" >> results.txt

validate_pdf.sh -c 7 -k 50 ../sample_classes/sctest.vec >> results.txt
echo "CASE5d" >> results.txt

validate_pdf.sh -n -c 7 -k 50 sq_reduced.vec >> results.txt
echo "CASE6a" >> results.txt

validate_pdf.sh -c 6 -k 500 -W 50 sq_reduced.vec >> results.txt
echo "CASE6b" >> results.txt

#browse_cluster_tree ../Landsat/cccd1.vec dum.cls < log.txt >> results.txt
#cls_comp_stats ../Landsat/cccd1.cls dum >> results.txt
#echo "CASE7" >> results.txt

