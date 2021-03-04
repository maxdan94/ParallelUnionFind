# Parallel Union Find
Feel free to use these lines as you wish. 
Parallel implementation of Union Find following that paper: https://papers-gamma.link/paper/193

## To compile

- "gcc Rem.c -o Rem -O3"
- "gcc RemPAR.c -o RemPAR -fopenmp -O3"

## To execute

- "./Rem edgelist.txt"

- "./RemPAR nthreads edgelist.txt"

"edgelist.txt" should contain one edge on each line "u v", where u and v are node IDs (unsigned int). 

Will print the number of edges in a spanning forest

## Initial contributors

Maximilien Danisch  
Mars 2021  
http://bit.ly/danisch  
maximilien.danisch@gmail.com


