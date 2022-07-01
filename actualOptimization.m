clc;
clear all;

len0=130;
A=[];
b=[];
Aeq = [];
beq = [];
lb = 100;
ub = 450;
nonlcon=@cCeqReturnFunction;
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
[le,fval,exitflag,output]=fmincon(@optimizationObjective,len0,A,b,Aeq,beq,lb,ub,nonlcon,options)