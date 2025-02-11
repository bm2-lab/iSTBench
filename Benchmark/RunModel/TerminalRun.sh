###Banksy
nohup Rscript Benchmark/RunModel/Run_Banksy.R \
-i /NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub/Data/BaristaSeq/sample_all_data/Slices_combind_data.RDS \
-o /NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub/Data/BaristaSeq/IntergrationRe \
-s Banksy -c 7 -d RDS -k 18 -l 0.2 \
> Data/BaristaSeq/IntergrationRe/Banksy.output &    

###CellCharter
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file /NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub/Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file /NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub/Data/BaristaSeq/IntergrationRe \
--sample CellCharter --nclust 7 --hvgs 2000 --runNormalization True \
--n_latent 15 --nhood_layers 2 \
> Data/BaristaSeq/IntergrationRe/CellCharter.output &

###CN
nohup python Benchmark/RunModel/Run_CN.py \
--input_file /NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub/Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file /NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub/Data/BaristaSeq/IntergrationRe \
--sample CN --nclust 7 \
> Data/BaristaSeq/IntergrationRe/CN.output &

###GraphST
nohup python Benchmark/RunModel/GraphST/Run_GraphST.py \
--input_file /NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub/Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file /NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub/Data/BaristaSeq/IntergrationRe \
--sample GraphST --nclust 7 --device cuda \
> Data/BaristaSeq/IntergrationRe/GraphST.output  &

###GraphST-PASTE
nohup python Benchmark/RunModel/Run_GraphST-PASTE.py \
--input_file /NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub/Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file /NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub/Data/BaristaSeq/IntergrationRe \
--sample GraphSTwithPASTE --nclust 7 --device cuda \
> Data/BaristaSeq/IntergrationRe/GraphSTwithPASTE.output  &

#MENDER
nohup python Benchmark/RunModel/Run_MENDER.py \
--input_file /NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub/Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file /NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub/Data/BaristaSeq/IntergrationRe \
--sample MENDER --nclust 7 --tech BaristaSeq \
> Data/BaristaSeq/IntergrationRe/MENDER.output &


#NicheCompass 
nohup python Benchmark/RunModel/Run_NicheCompass.py \
--input_file /NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub/Data/BaristaSeq/sample_data \
--output_file /NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub/Data/BaristaSeq/IntergrationRe \
--sample NicheCompass --nclust 7 \
> Data/BaristaSeq/IntergrationRe/NicheCompass.output &

#Spado
nohup Rscript Benchmark/RunModel/Run_Spado.R \
-i /NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub/Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
-o /NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub/Data/BaristaSeq/IntergrationRe \
-s Spado -v 2000 -c 7 -n TRUE \
> Data/BaristaSeq/IntergrationRe/Spado.output &    

