# Generalized Damped Newton Method for Solving Lasso - GDNM

### In [[kmpt21]](https://arxiv.org/abs/2101.10555), we propose a new algorithm for solving Lasso regession which is the optimization problem of the form
![image](https://user-images.githubusercontent.com/69850027/107838172-8fba9080-6d72-11eb-8126-0d2d833331b7.png)
![image](https://user-images.githubusercontent.com/69850027/107067479-55b72080-67ad-11eb-937e-1a5d55cee152.png)
![image](https://user-images.githubusercontent.com/69850027/107069564-23f38900-67b0-11eb-8e76-65f5ac3084ee.png)
### The Matlab code for this algorithm can be found in the file GDNM.m or below.
### Input 
```rb
function x=lasso_GDNM(A,b,mu)
```
### For simplicity, we choose \beta=1/2 and \sigma=0.1. The choice of y^0 can be found in the iterating step below. 

### Find gamma.
```rb
    ATA=A'*A;
    gamma=0.5/eigs(ATA,1);
```
## Calculate Proximal mapping
```rb
    function [value,zeros_index]=prox_of_function(gamma,y,mu)
        value=[];
        zeros_index=[];
        for i=1:length(y)
            if y(i)>mu*gamma
                value=[value;y(i)-mu*gamma];
            elseif y(i)<-mu*gamma
                value=[value;y(i)+mu*gamma];
            else 
                value=[value;0];
                zeros_index=[zeros_index;i];
            end
        end
    end   
 ```
 ### Calculate matrices Q, P, and vector c
```rb
    Q=inv(eye(size(ATA))-gamma*ATA);
    P=Q-eye(size(Q));
    c=-gamma*Q*transpose(A)*b;
```
### Calculate value of Lasso function, value of function phi (1) and its gradient (2)
```rb    
    function value=value_of_function(A,b,mu,x)
        value=0.5*norm(A*x-b)^2+mu*norm(x,1);
    end
    function value_of_phi=phi(A,b,gamma,mu,y)
        [prox_y,zeros_index]=prox_of_function(gamma,y,mu);
        value_of_phi=0.5*transpose(P*y)*y+transpose(c)*y+gamma*mu*norm(prox_y,1)+0.5*norm(y-prox_y)^2;
    end
    function grad=finding_gradient(A,b,gamma,mu,y)
        [prox_y,zeros_index]=prox_of_function(gamma,y,mu);
        grad=Q*y-prox_y+c;
    end
```
### Find direction dk, step size tk
```rb
    function dk=finding_dk(A,b,gamma,mu,y)
        grad=finding_gradient(A,b,gamma,mu,y);
        [prox_y,zeros_index]=prox_of_function(gamma,y,mu);
        X=P;
        for i=zeros_index
            X(i,:)=Q(i,:);
        end
        dk=-linsolve(X,grad);
    end
    function tau=finding_tk(A,b,gamma,mu,y)
        tau=1;
        dk=finding_dk(A,b,gamma,mu,y);
        grad=finding_gradient(A,b,gamma,mu,y);
        while phi(A,b,gamma,mu,y+tau*dk)-phi(A,b,gamma,mu,y)-0.1*tau*grad'*dk>0
            tau=tau/2;
        end
    end
 ```
 ### Start iterating
 ```rb
    y=zeros(length(c),1);
    iter=0;
    while iter<100
        y=y+finding_tk(A,b,gamma,mu,y)*finding_dk(A,b,gamma,mu,y);
        iter=iter+1
        value=0.5*norm(A*(Q*y+c)-b)^2+mu*norm(Q*y+c,1)
    end
 ```
### Translate y to x
 ```rb
    x=Q*y+c;
end
```
### The entire code can be found in file GDNM.m

### Let us do a simple example from [[Example 7.5,kmp20]](https://arxiv.org/abs/2009.10551) for illustration. In this example, we choose 
![image](https://user-images.githubusercontent.com/69850027/107070744-b6e0f300-67b1-11eb-8ab7-5f10f9a7eaa0.png)
### The implementations for this are as follows
 ```rb
A=[4/7,3,6;12/7,2,-3;6/7,-6,2;0,24,0];
b=[104/49;347/49;-649/49;0];
mu=1/3;
x=lasso_GDNM(A,b,mu)
```
### Reader can find implementations for this example via Octave Online in the link [GDNM Example](http://bit.ly/GDMN_example).
## References
### [kmp20] Khanh, P. D., Mordukhovich, B. S., Phat, V. T.: A Generalized Newton Method for Subgradient Systems, submitted, https://arxiv.org/abs/2009.10551
### [kmpt21] Khanh, P. D., Mordukhovich, B. S., Phat, V. T., Tran, D. B.: Generalized Damped Newton Algorithms in Nonsmooth Optimization with Applications to Lasso Problems, submitted, https://arxiv.org/abs/2101.10555 
