for s in 020 035 050; do
  for m in 050 100 150 200; do
    name=dg_s${s}m${m}
    echo " "
    echo " ----------------------------------------------------  "
    echo ${name}
    python generate_mesh.py -n ${name} -g doughnut -size ${s} -margin ${m} -highr 3 -lowr 20 -p || echo "Error or done."
    echo " "
    echo " ----------------------------------------------------  "
    echo " "
  done
done