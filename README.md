# FAQ

* Canonical reads: In a sequence read, a read cannot be said to have originated from a particular strand of DNA. So there are two ways to resolve the choice between a read and its complement. Discussion on this [forum](https://www.biostars.org/p/153170/) that method 2 should be used.
    1. Take the lexicographically smallest of read and canonical form
    2. Consider both cases as unique
    
 Right now , we are not considering them unique. I'll check out what perga does . 
 
* PERGA CODE: https://github.com/zhuxiao/PERGA
 
* In the synchronization step for binning, While merging the bins, we need to know which bins will stay on which node. For that , we need to know the total no. of bins and the total no of nodes. To know the total no. of bins , we need another global communication. Is there any other way to do this?

* Data sets : https://www.ebi.ac.uk/ena

 
