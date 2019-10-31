#! /bin/bash

#./buildmulti.sh hkdr.valor/listSK_stat     &> /dev/null
#./buildmulti.sh hkdr.valor/listT2HK_stat   &> /dev/null
#./buildmulti.sh hkdr.valor/listSKT2HK_stat &> /dev/null
#
#./buildmulti.sh hkdr.valor/listSK_syst     &> /dev/null
#./buildmulti.sh hkdr.valor/listT2HK_syst   &> /dev/null
#./buildmulti.sh hkdr.valor/listSKT2HK_syst &> /dev/null
#
#./buildmulti.sh asimov.valor/listSK_stat     &> /dev/null
#./buildmulti.sh asimov.valor/listT2HK_stat   &> /dev/null
#./buildmulti.sh asimov.valor/listSKT2HK_stat &> /dev/null
#
#./buildmulti.sh asimov.valor/listSK_syst     &> /dev/null
#./buildmulti.sh asimov.valor/listT2HK_syst   &> /dev/null
#./buildmulti.sh asimov.valor/listSKT2HK_syst &> /dev/null


###only test points

./buildmulti.sh hkdr_nh_ih.valor/testSK_stat     
./buildmulti.sh hkdr_nh_ih.valor/testT2HK_stat   
./buildmulti.sh hkdr_nh_ih.valor/testSKT2HK_stat 

./buildmulti.sh hkdr_ih_nh.valor/testSK_stat     
./buildmulti.sh hkdr_ih_nh.valor/testT2HK_stat   
./buildmulti.sh hkdr_ih_nh.valor/testSKT2HK_stat 

./buildmulti.sh asimov_nh_ih.valor/testSK_stat     
./buildmulti.sh asimov_nh_ih.valor/testT2HK_stat   
./buildmulti.sh asimov_nh_ih.valor/testSKT2HK_stat 

./buildmulti.sh asimov_ih_nh.valor/testSK_stat     
./buildmulti.sh asimov_ih_nh.valor/testT2HK_stat   
./buildmulti.sh asimov_ih_nh.valor/testSKT2HK_stat 
