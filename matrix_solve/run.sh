set -x

bsub -b -I -q q_cpc -n 6 -np 1 -cgsp 64 -host_stack 1024 -share_size 15000 -cache_size 128  ./matrix_solve
