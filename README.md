#  Lifted Disjoint Paths with Application in Multiple Object Tracking

[![PWC](https://img.shields.io/endpoint.svg?url=https://paperswithcode.com/badge/lifted-disjoint-paths-with-application-in-1/multi-object-tracking-on-2d-mot-2015)](https://paperswithcode.com/sota/multi-object-tracking-on-2d-mot-2015?p=lifted-disjoint-paths-with-application-in-1) [![PWC](https://img.shields.io/endpoint.svg?url=https://paperswithcode.com/badge/lifted-disjoint-paths-with-application-in-1/multi-object-tracking-on-mot16)](https://paperswithcode.com/sota/multi-object-tracking-on-mot16?p=lifted-disjoint-paths-with-application-in-1)  	
[![PWC](https://img.shields.io/endpoint.svg?url=https://paperswithcode.com/badge/lifted-disjoint-paths-with-application-in-1/multi-object-tracking-on-mot17)](https://paperswithcode.com/sota/multi-object-tracking-on-mot17?p=lifted-disjoint-paths-with-application-in-1)

This is the official implementation of our **ICML 2020** paper *Lifted Disjoint Paths with Application in Multiple Object Tracking* ([Andrea Hornakova](https://www.mpi-inf.mpg.de/departments/computer-vision-and-machine-learning/people/andrea-hornakova), [Roberto Henschel](http://www.tnt.uni-hannover.de/staff/henschel/), [Bodo Rosenhahn](http://www.tnt.uni-hannover.de/en/staff/rosenhahn/), [Paul Swoboda](https://www.mpi-inf.mpg.de/departments/computer-vision-and-machine-learning/people/paul-swoboda/)) [https://arxiv.org/abs/2006.14550].


![Tracking Result](data/output.gif)

We provide the solver implemented in C++, along with a Python wrapper.
The tracker is implemented in Matlab.



## Evaluation
Using the features explained in our paper, we achieve the following results on MOT17:


|           | MOTA         | IDF1           |       FP     |     FN     |     IDs      |
|  :---:    | :---:        |     :---:      |    :---:     | :---:      |    :---:     |
| **Train** |     67.0     |     72.4       |    2655      |   107803   |     791      |
| **Test**  |     60.5     |     65.6       |    14966     | 206619     |     1189     |


Note that all results on the training set have been calculated in a leave-one-out fashion so that values are actually meaningful. 

## Citation
If you use our work in your research, please cite our publication:
```
    @InProceedings{lifted_disjoint_paths_2020_ICML,
    author={Andrea Hornakova and Roberto Henschel and Bodo Rosenhahn and Paul Swoboda},
    title={Lifted Disjoint Paths with Application in Multiple Object Tracking},
    booktitle = {The 37th International Conference on Machine Learning (ICML)},
    month = {July},
    year = {2020}
}
```



