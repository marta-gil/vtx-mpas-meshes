for s in 010 020 030; do
  for m in 050 100 150; do
    name=dg1s${s}m${m}
    echo " "
    echo " ----------------------------------------------------  "
    echo ${name}
    python generate_mesh.py -n ${name} -g doughnut -size ${s} -margin ${m} -highr 3 -lowr 20 -p || echo "Error or done."
    echo " "
    echo " ----------------------------------------------------  "
    echo " "
  done
done