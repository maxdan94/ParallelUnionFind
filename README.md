# Parallel Union Find
Feel free to use these lines as you wish. 
Parallel implementation of Union Find following that paper: https://papers-gamma.link/paper/193

## To compile

- "gcc Rem.c -o Rem -O3"
- "gcc Rem_lock.c -o Rem_lock -fopenmp -O3"
- "gcc Rem_verif.c -o Rem_verif -fopenmp -O3"



## To execute

- "./Rem edgelist.txt"
- "./Rem_lock nthreads edgelist.txt"
- "./Rem_verif nthreads edgelist.txt"

"edgelist.txt" should contain one edge on each line "u v", where u and v are node IDs (unsigned int). 

Will print the number of edges in a spanning forest (Rem_verif can output wrong result: a few more edges than what should be)

## Initial contributors

Maximilien Danisch  
Mars 2021  
http://bit.ly/danisch  
maximilien.danisch@gmail.com


