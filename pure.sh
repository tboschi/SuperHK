#! /bin/bash

#./prepare_systematics.sh -r /home/tboschi/data/errorstudy/0   DataFit_Postfit_2018_final_v1_180907_sk_eb_valor_order_all_plus_scc.root SKJointErrorMatrix2018_Total_fitqun_v4_16thAug2017_VALOR_order.root
#./prepare_systematics.sh -r /home/tboschi/data/errorstudy/1a  banff_00_nonuenumucc_remove.root SKJointErrorMatrix2018_Total_fitqun_v4_16thAug2017_VALOR_order.root 
#./prepare_systematics.sh -r /home/tboschi/data/errorstudy/1b  banff_01_nonuenumucc_remove.root SKJointErrorMatrix2018_Total_fitqun_v4_16thAug2017_VALOR_order.root 
#./prepare_systematics.sh -r /home/tboschi/data/errorstudy/2a  banff_02_nonuebarnumubarcc_remove.root SKJointErrorMatrix2018_Total_fitqun_v4_16thAug2017_VALOR_order.root
#./prepare_systematics.sh -r /home/tboschi/data/errorstudy/2b  banff_03_nonuebarnumubarcc_remove.root SKJointErrorMatrix2018_Total_fitqun_v4_16thAug2017_VALOR_order.root
#./prepare_systematics.sh -r /home/tboschi/data/errorstudy/6a  banff_nonu2p2h_remove.root SKJointErrorMatrix2018_Total_fitqun_v4_16thAug2017_VALOR_order.root
#./prepare_systematics.sh -r /home/tboschi/data/errorstudy/7a  banff_nonubar2p2h_remove.root SKJointErrorMatrix2018_Total_fitqun_v4_16thAug2017_VALOR_order.root
#./prepare_systematics.sh -r /home/tboschi/data/errorstudy/67  banff_nonu2p2h_nonubar2p2h_remove.root SKJointErrorMatrix2018_Total_fitqun_v4_16thAug2017_VALOR_order.root
#./prepare_systematics.sh -r /home/tboschi/data/errorstudy/8   banff_increasednueflux.root SKJointErrorMatrix2018_Total_fitqun_v4_16thAug2017_VALOR_order.root
#./prepare_systematics.sh -r /home/tboschi/data/errorstudy/9   banff_06.root SKJointErrorMatrix2018_Total_fitqun_v4_16thAug2017_VALOR_order.root
#./prepare_systematics.sh -r /home/tboschi/data/errorstudy/10  banff_07.root SKJointErrorMatrix2018_Total_fitqun_v4_16thAug2017_VALOR_order.root
#./prepare_systematics.sh -r /home/tboschi/data/errorstudy/11a DataFit_Postfit_2018_final_v1_180907_sk_eb_valor_order_all_plus_scc.root skdetfsi_escale29.root
#./prepare_systematics.sh -r /home/tboschi/data/errorstudy/11b DataFit_Postfit_2018_final_v1_180907_sk_eb_valor_order_all_plus_scc.root skdetfsi_escale19.root
#./prepare_systematics.sh -r /home/tboschi/data/errorstudy/NC  banff_14.root SKJointErrorMatrix2018_Total_fitqun_v4_16thAug2017_VALOR_order.root
#./prepare_systematics.sh -r /home/tboschi/data/errorstudy/flux_lukas  banff_extraflux.root SKJointErrorMatrix2018_Total_fitqun_v4_16thAug2017_VALOR_order.root
./prepare_systematics.sh -r /home/tboschi/data/errorstudy/nuenorm1_corr norm_nue1_nuebar1_corr.root
./prepare_systematics.sh -r /home/tboschi/data/errorstudy/nuenorm2_corr norm_nue2_nuebar2_corr.root
./prepare_systematics.sh -r /home/tboschi/data/errorstudy/nuenorm3_corr norm_nue3_nuebar3_corr.root
./prepare_systematics.sh -r /home/tboschi/data/errorstudy/nuenorm4_corr norm_nue4_nuebar4_corr.root
./prepare_systematics.sh -r /home/tboschi/data/errorstudy/nuenorm5_corr norm_nue5_nuebar5_corr.root

./prepare_systematics.sh -r /home/tboschi/data/errorstudy/nuenorm1_anti norm_nue1_nuebar1_anticorr.root
./prepare_systematics.sh -r /home/tboschi/data/errorstudy/nuenorm2_anti norm_nue2_nuebar2_anticorr.root
./prepare_systematics.sh -r /home/tboschi/data/errorstudy/nuenorm3_anti norm_nue3_nuebar3_anticorr.root
./prepare_systematics.sh -r /home/tboschi/data/errorstudy/nuenorm4_anti norm_nue4_nuebar4_anticorr.root
./prepare_systematics.sh -r /home/tboschi/data/errorstudy/nuenorm5_anti norm_nue5_nuebar5_anticorr.root
