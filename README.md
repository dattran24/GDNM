# GDNM
## Generalized Damped Newton Method for solving Lasso

### In this document, we propose a new algorithm for solving Lasso regession, i. e. the optimization problem of the form

![h1](https://user-images.githubusercontent.com/69850027/106015925-f52d3280-608c-11eb-8c10-de908dba6c94.png)

![h2](https://user-images.githubusercontent.com/69850027/106018853-fb70de00-608f-11eb-9612-42b746aac2d2.png)


### This algorithm, as proved by numerical experiments in https://arxiv.org/abs/2101.10555, outperformed well-known algorithms such as FISTA, APG, SSNAL, ADMM in the case m>=n. To proceed, we need the notation of proximal mapping, as suggested in the paper above.

![h3](https://user-images.githubusercontent.com/69850027/106018740-ded4a600-608f-11eb-9343-8f1a44a5487c.png)

![h4](https://user-images.githubusercontent.com/69850027/106019826-142dc380-6091-11eb-9c64-55764e3bf446.png)

