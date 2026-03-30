echo "N,R,A,C,w,g" > results.csv
~/git/monolayergrowth/results/postprocessing/metrics cells.csv boundary.csv neighbors.csv >> results.csv

#Input:
#cells.csv : CSV file containing the x,y,g,n table from your simulation (x,y are centroid coordinates per cell in length unit "R" (equilibrium cell radius) and not pixels if you simulate on-lattice, g is 0 for inhibited and 1 for growing, n is the number of neighbors).

#Output:

#boundary.csv : Optional CSV file containing the ordered x,y table of all boundary cells
#neighbors.csv : Optional CSV file containing the neighbor number distribution
#results.csv : CSV list of total number of cells (N), total tissue radius (R), total tissue area (A), tissue circumference (C), boundary width (w), fraction of growing cells (g)

#To compute the boundary roughness, calculate C/(2*pi*sqrt(A/pi)) or w/R (we have two ways of quantifying it).
#To batch process a series of files (for example a time series), run something like:
#echo "N,R,A,C,w,g" > results.csv
#for i in {0..1000}
#do
#    ./metrics cells$i.csv boundary$i.csv neighbors$i.csv >> results.csv
#done

#This will add one line per timepoint to results.csv.

