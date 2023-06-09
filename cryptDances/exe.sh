avg_time() {
    #
    # usage: avg_time n command ...
    #
    n=$1; shift
    (($# > 0)) || return                   # bail if no command given
    for ((i = 0; i < n; i++)); do
        { time -p "$@" &>/dev/null; } 2>&1 # ignore the output of the command
                                           # but collect time's output in stdout
    done | awk '
        /real/ { real = real + $2; nr++ }
        /user/ { user = user + $2; nu++ }
        /sys/  { sys  = sys  + $2; ns++}
        END    {
                 if (nr>0) printf("real %f\n", real/nr);
                 if (nu>0) printf("user %f\n", user/nu);
                 if (ns>0) printf("sys %f\n",  sys/ns)
               }'
}
make
#avg_time 5 mpirun --mca btl_tcp_if_include eno1 -np 32 -host 10.191.12.53:8,10.191.12.54:8,10.191.12.55:8,172.31.102.56:8 ./crypt_dances_explorer
avg_time 5 mpirun --mca btl_tcp_if_include eno1 -np 8 -host 10.191.12.54:8 ./crypt_dances_explorer