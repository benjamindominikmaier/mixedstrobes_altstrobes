#!/bin/bash


genome1=$1
genome2=$2

IFS=$'\n'       # make newlines the only separator
# set -f          # disable globbing


#######################################################
#######################################################
########## INFORMATION ABOUT THE SCRIPT ###############
#######################################################
#######################################################

# RUN scripts e.g. as:   ./genome_mappings_eval.sh <genome1> <genome2>

# The mummerplot lines can be activated for plotting dotplots 

# WARNING: The genome_mapping_metrics.py script should not be used 
# for large outfiles (e.g., files with more than 10M-20M matches as
# it uses too much memory and becomes slow. 

for k in 30 100 500
do
	# mummer MEM
	/usr/bin/time -f '%e real, %M  maximum resident set size'  mummer -F -maxmatch -l $k -b  $genome1 $genome2 > tmp.tsv 2> runtime.txt
	echo -n "MUMmer & MEM & " $k " & " 
	python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv

	# # mummer MUM
	/usr/bin/time -f '%e real, %M  maximum resident set size'  mummer -F  -l $k -mum -b $genome1 $genome2 > tmp.tsv 2> runtime.txt
	echo -n "MUMmer & MUM & " $k " & " 
	python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

done

echo "Start StrobeMap"

# StrobeMap
/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -k 30 -v 31 -c kmers -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap  & 30 & 100 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 15 -v 1 -w 100 -c randstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap  & (2, 15, 1, 100) & 100 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 15 -v 16 -w 100 -c randstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap  & (2, 15, 16, 100) & 100 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv

#/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 3 -k 10 -v 11 -w 100 -c randstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
#echo -n "./StrobeMap & (3, 10, 11, 100) & 100 & " 
#python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv



/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 30 -v 31 -w 100 -c randstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt  1> stdout.txt
echo -n "./StrobeMap & (2, 30, 31, 100) & 100 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 3 -k 30 -v 31 -w 100 -c randstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt  1> stdout.txt
echo -n "./StrobeMap & (3, 30, 31, 100) & 100 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 30 -v 200 -w 400 -c randstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap & (2, 30, 200, 400) & 100 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 3 -k 30 -v 200 -w 400 -c randstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap & (3, 30, 200, 400) & 100 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 30 -v 500 -w 600 -c randstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap & (2, 30, 500, 600) & 100 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


#######################################################
#######################################################
#######################################################


# StrobeMap
/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 15 -v 1 -f 30 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap  & (2, 15, 1, 100) & 30 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 15 -v 16 -f 30 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap  & (2, 15, 16, 100) & 30 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv

# /usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 3 -k 10 -v 11 -f 30 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
# echo -n "./StrobeMap & (3, 10, 11, 100) & 30 & " 
# python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv



/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 30 -v 31 -f 30 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt  1> stdout.txt
echo -n "./StrobeMap & (2, 30, 31, 100) & 30 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 3 -k 30 -v 31 -f 30 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt  1> stdout.txt
echo -n "./StrobeMap & (3, 30, 31, 100) & 30 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 30 -v 200 -f 30 -w 400 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap & (2, 30, 200, 400) & 30 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 3 -k 30 -v 200 -f 30 -w 400 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap & (3, 30, 200, 400) & 30 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 30 -v 500 -f 30 -w 600 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap & (2, 30, 500, 600) & 30 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


#######################################################
#######################################################
#######################################################


# StrobeMap
/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 15 -v 1 -f 45 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap  & (2, 15, 1, 100) & 45 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 15 -v 16 -f 45 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap  & (2, 15, 16, 100) & 45 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv

# /usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 3 -k 10 -v 11 -f 45 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
# echo -n "./StrobeMap & (3, 10, 11, 100) & 45 & " 
# python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv



/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 30 -v 31 -f 45 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt  1> stdout.txt
echo -n "./StrobeMap & (2, 30, 31, 100) & 45 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 3 -k 30 -v 31 -f 45 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt  1> stdout.txt
echo -n "./StrobeMap & (3, 30, 31, 100) & 45 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 30 -v 200 -f 45 -w 400 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap & (2, 30, 200, 400) & 45 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 3 -k 30 -v 200 -f 45 -w 400 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap & (3, 30, 200, 400) & 45 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 30 -v 500 -f 45 -w 600 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap & (2, 30, 500, 600) & 45 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


#######################################################
#######################################################
#######################################################


# StrobeMap
/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 15 -v 1 -f 60 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap  & (2, 15, 1, 100) & 60 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 15 -v 16 -f 60 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap  & (2, 15, 16, 100) & 60 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv

# /usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 3 -k 10 -v 11 -f 60 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
# echo -n "./StrobeMap & (3, 10, 11, 100) & 60 & " 
# python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv



/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 30 -v 31 -f 60 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt  1> stdout.txt
echo -n "./StrobeMap & (2, 30, 31, 100) & 60 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 3 -k 30 -v 31 -f 60 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt  1> stdout.txt
echo -n "./StrobeMap & (3, 30, 31, 100) & 60 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 30 -v 200 -f 60 -w 400 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap & (2, 30, 200, 400) & 60 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 3 -k 30 -v 200 -f 60 -w 400 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap & (3, 30, 200, 400) & 60 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 30 -v 500 -f 60 -w 600 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap & (2, 30, 500, 600) & 60 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


#######################################################
#######################################################
#######################################################


# StrobeMap
/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 15 -v 1 -f 60 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap  & (2, 15, 1, 100) & 75 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 15 -v 16 -f 60 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap  & (2, 15, 16, 100) & 75 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv

# /usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 3 -k 10 -v 11 -f 60 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
# echo -n "./StrobeMap & (3, 10, 11, 100) & 75 & " 
# python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv



/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 30 -v 31 -f 75 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt  1> stdout.txt
echo -n "./StrobeMap & (2, 30, 31, 100) & 75 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 3 -k 30 -v 31 -f 75 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt  1> stdout.txt
echo -n "./StrobeMap & (3, 30, 31, 100) & 75 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 30 -v 200 -f 75 -w 400 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap & (2, 30, 200, 400) & 75 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 3 -k 30 -v 200 -f 75 -w 400 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap & (3, 30, 200, 400) & 75 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 30 -v 500 -f 75 -w 600 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap & (2, 30, 500, 600) & 75 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


#######################################################
#######################################################
#######################################################

# StrobeMap
/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 15 -v 1 -f 90 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap  & (2, 15, 1, 100) & 90 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 15 -v 16 -f 90 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap  & (2, 15, 16, 100) & 90 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv

# /usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 3 -k 10 -v 11 -f 90 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
# echo -n "./StrobeMap & (3, 10, 11, 100) & 90 & " 
# python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2 --collinear_matches_out coll_sol.tsv



/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 30 -v 31 -f 90 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt  1> stdout.txt
echo -n "./StrobeMap & (2, 30, 31, 100) & 90 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 3 -k 30 -v 31 -f 90 -w 100 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt  1> stdout.txt
echo -n "./StrobeMap & (3, 30, 31, 100) & 90 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 30 -v 200 -f 90 -w 400 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap & (2, 30, 200, 400) & 90 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 3 -k 30 -v 200 -f 90 -w 400 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap & (3, 30, 200, 400) & 90 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

/usr/bin/time -f '%e real, %M  maximum resident set size'  ./StrobeMap -n 2 -k 30 -v 500 -f 90 -w 600 -c mixedrandstrobes -S -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "./StrobeMap & (2, 30, 500, 600) & 90 & " 
python3 genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


#######################################################
#######################################################
#######################################################

