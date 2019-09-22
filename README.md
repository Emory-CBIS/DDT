![alt text](gui/cbis.png "CBIS")

# Difference Degree Test toolbox 

The Difference Degree Test (DDT) is a two stage method to detect regions incident to a statistically significant number of differentially weighted edges (DWEs). In the phase, we select a data-adaptive threshold to identify the DWEs followed by a statistical test for the number of DWEs incident to each brain region. The key to our procedure the Hirscheberger-Qi-Steuer (Hirschberger et al., 2007) algorithm, which is a computationally efficient algorithm for generating random null networks that replicate statistical properties of the observed difference network. 

Further details on the DDT, please see our paper (Higgins et al., 2019). 



## Usage

Store the correlation matrices in a cell array.  All transformations will be performed by DDT.  The covariate matrix must also be a matrix with the group membership variable as the first column.  For categorical variables with k>2 levels, please create k-1 dummy variables.

```
function [D,pvals,vec,numDWE] = DDT(corrmats,covariate,adjust,method_null,M)
% INPUT
%   corrmats:     an array containing correlation matrices for all subjects
%   covariate:    design matrix (first column must be the group membership)
%   adjust:       binary indicator to adjust (=1) or not adjust (=0) for
%                 variables in covariate
%   method_null:  'eDDT' (empirical threshold) or 'aDDT' (theoretical
%                 threshold
%   M:            number of HQS nulls to generate

%OUTPUT
%   D:            Difference network (ignore diagonal)
%   pvals:        p value for the test at each node
%   vec:          regions incident to significant num of DWEs
%   numDWE:       number of DWEs incident to each node


```



## Version

DDT is currently in version 0.0.1

## License

This toolbox is licensed under the MIT License - see the LICENSE file for details

## Contact

Please send comments and bug reports to: ihiggin@emory.edu

## References

Higgins, I.A., Kundu, S., Choi, K.S., Mayberg, H.S. and Guo, Y. (2019). A difference degree test for comparing brain networks. Human brain mapping, 40(15):4518-4536.

Hirschberger, M., Qi, Y., and Steuer, R. E. (2007). Randomly generating portfolio-selection covariance matrices with specified distributional characteristics. European Journal of Operational Research, 177(3):1610–1625.
