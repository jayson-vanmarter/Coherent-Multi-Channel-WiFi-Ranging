# Coherent-Multi-Channel-WiFi-Ranging

This repository contains the data and code related to our paper "Coherent Multi-Channel WiFi Ranging for Next Generation Positioning." This work demonstrates the potential of multi-channel WiFi ranging using timestamp and CFR measurements supporting Next Generation Positioning. Utilizing our proposed methods and collected data in real-world indoor evironments, we demonstrate a 90th percentile ranging error of 9 cm in line-of-sight (LOS) conditions for stitched bandwidths of at least 320 MHz.

Data set:  
https://utdallas.box.com/v/WiFi-Ranging-Testbed-Data

If you use our code or data in your research, please cite our paper:  
```
@article{VanMarter-WiFiRanging,
  author = {Van Marter, Jayson P. and Ben-Shachar, Matan and Alpert, Yaron and Dabak, Anand G. and Al-Dhahir, Naofal and Torlak, Murat},
  journal = {}, 
  title = {Coherent Multi-Channel WiFi Ranging for Next Generation Positioning}, 
  year = {2023},
}
```

## Quick Start

1. Download the provided data set
2. Set the trial names and distances in the _process_data_trials.m_.
3. Set the packet type and bandwidth in _process_data_trials.m_.
4. Run _process_data_trials.m_.

An example plot across multiple bandwidths using results from running _process_data_trials.m_ can be created by running _plot_multiple_results.m_.

# License
[GPL 3.0](https://choosealicense.com/licenses/gpl-3.0/)
