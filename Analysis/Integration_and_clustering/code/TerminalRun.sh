#1.1 Integration_metric.py
#BaristaSeq
nohup python Intergration_metric.py \
--input_file /home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark/BaristaSeq/IntergrationRe \
--input_path /home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_path /home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark/BaristaSeq/IntergrationRe/Metric/Intergration_value.csv \
> /home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark/BaristaSeq/IntergrationRe/Metric/Intergration_metric.output &

#1.2 Domain_plotting.R
#BaristaSeq
nohup /usr/lib/R/bin/Rscript Domain_plot.R \
/home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark/BaristaSeq/IntergrationRe/Metric/domain_plot.txt \
/home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark/BaristaSeq/IntergrationRe/Metric \
> /home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark/BaristaSeq/IntergrationRe/Metric/domain_plot.output&
