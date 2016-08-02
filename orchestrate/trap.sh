function finish {
    # re-start service
   ./CMD
}
while :
    do
        trap finish EXIT
        ./CMD
done

