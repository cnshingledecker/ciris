function finish {
    # re-start service
   ./RCMD
}
while :
    do
        trap finish EXIT
        ./RCMD
done

