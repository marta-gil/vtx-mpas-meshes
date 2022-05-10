names=`ls /home/marta/PycharmProjects/vtx-mpas-meshes/data`
names2=`echo $names | sed 's/ /,/g'`

python compare_meshes.py -gs ${names2} -o compared.pdf
