# Generalized Damped Newton Method for solving Lasso - GDNM

### In this document, we propose a new algorithm for solving Lasso regession, i. e. the optimization problem of the form

![image](https://user-images.githubusercontent.com/69850027/107066957-9e220e80-67ac-11eb-949e-1e166d81b13c.png)
![image](https://user-images.githubusercontent.com/69850027/107067479-55b72080-67ad-11eb-937e-1a5d55cee152.png)
![image](https://user-images.githubusercontent.com/69850027/107067681-a3338d80-67ad-11eb-8761-4a20e1b4aec7.png)
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

## Reference
### [kmpt21] Khanh, P. D., Mordukhovich, B. S., Phat, V. T., Tran, D. B.: Generalized Damped Newton Algorithms in Nonsmooth Optimization with Applications to Lasso Problems, submitted, https://arxiv.org/abs/2101.10555 
