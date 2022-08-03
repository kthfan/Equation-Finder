# Equation-Finder
Given $x,y$, find a equation $f$ such that $y = f(x)$.
This equation finding program is implemented in genetic algorithm. 

# Interface
<img width="750" src="https://github.com/kthfan/Equation-Finder/blob/main/screenshot/interface.jpg"></img>

# Usage
## Basic usage
Press "read variables" button to load variables from .mat file. The mat file must include $x,y$, where $x\in \mathbb{R}^{N\times C}$,  $x\in \mathbb{R}^{N\times 1}$, $N$ is number of samples and $C$ is number of categories.  
Press "solve" button to search the nearly approximate equation according to loaded variables $x,y$.

<img width="750" src="https://github.com/kthfan/Equation-Finder/blob/main/screenshot/solved.jpg"></img>

The "loss" edit field shows the mean squared error $\rVert y-\hat{f}(x)\lVert_2$, where $\hat{f}$ is estimated equation of $f$. 
The "equation" edit field display a matlab executable statement refers to $\hat{f}$.
A beautify estimated equation is showed on the bottom of UI.  

The "export model" button allow the trained model saved in the storage and "import model" button loads saved model into program.
The "reset" button reset the model and clear all variables that store training information.

# References
1. Neapolitan, R., & Naimipour, K. (2010). Foundations of algorithms 5/e. Jones & Bartlett Publishers. In chapter 10.3.
