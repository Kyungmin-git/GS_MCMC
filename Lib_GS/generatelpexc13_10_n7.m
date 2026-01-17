function [candidatet0,candidateqn_n,gastype,logJacobian,logproposal]=generatelpexc13_10_n7(currentt0,currentqn_n)
global tau gastype
Ncurrent=length(currentt0);
gastype=randsample(4,1);

if Ncurrent<101
gastype=randsample(3,1);
elseif Ncurrent>1999
    gastype=randsample([1,2,4],1);
else
    gastype=randsample(4,1);
end

Stepsize=0.02*abs(normrnd(0,1));
Stepsize2=1;
alpha=1;
logJacobian=0;
logproposal=0;
if  gastype==1   
    Nd=max(1,round(Stepsize*size(currentt0,2)));
    ixshift=randperm(size(currentt0,2),Nd);
    currentt0(ixshift)=currentt0(ixshift)+normrnd(0,0.04,[1,Nd]);
    currentt0=abs(currentt0);
    idxt0=find(currentt0>tau);
    currentt0(idxt0)=2*tau-currentt0(idxt0);
    candidatet0=currentt0;
    candidateqn_n=currentqn_n;
    gas_matrix=[candidatet0;candidateqn_n]; 
    gas_matrix_vert=gas_matrix';
    gas_matrix_vert=sortrows(gas_matrix_vert);
    candidatet0=gas_matrix_vert(:,1)';
    candidateqn_n=gas_matrix_vert(:,2)';
    
    
   
elseif gastype ==2
    Nd=max(1,round(Stepsize*size(currentt0,2)));
    qxshift=randperm(size(currentt0,2),Nd);
    candidateqn_n=currentqn_n;
    for i=1:Nd
    candidateqn_n(qxshift(i))=currentqn_n(qxshift(i))+normrnd(0,0.00015,1);
    end
    candidateqn_n=abs(candidateqn_n);
    candidateqn_n=candidateqn_n/sum(candidateqn_n);
    qn_dif=sum(candidateqn_n(qxshift)-currentqn_n(qxshift));
    nj=length(currentqn_n)-Nd;
    Jacobian=(1-qn_dif)^(nj);
    logJacobian=log(Jacobian);
    %{
    re_candidateqn=candidateqn;
    resid_qn=sum(candidateqn)-sum(currentqn);
    re_candidateqn(qxshift)=0;
    sum_re_candidateqn=sum(re_candidateqn);
    re_candidateqn=re_candidateqn*(sum_re_candidateqn+resid_qn)/sum_re_candidateqn;
    re_candidateqn(qxshift)=candidateqn(qxshift);
    candidateqn=re_candidateqn;
    logcandidateqn=log10(candidateqn);
    %}
    candidatet0=currentt0;
    
elseif gastype ==3 %% birth
    if Ncurrent==1999
        birthn=1;
    else
    birthn=randperm(2,1);
    end
    Nt=length(currentt0);
    t0_new=rand(1,birthn)*tau;
    candidatet0=[currentt0 t0_new]; 
    qn_n_new=zeros(1,birthn);
    for i=1:birthn
       new_x=betarnd(alpha,(Nt+birthn-1)*alpha);
       qn_n_new(i)=new_x;
    end
    currentqn_n=currentqn_n*(1-sum(qn_n_new));
    candidateqn_n=[currentqn_n qn_n_new];
    gas_matrix=[candidatet0;candidateqn_n]; 
    gas_matrix_vert=gas_matrix';
    gas_matrix_vert=sortrows(gas_matrix_vert);
    candidatet0=gas_matrix_vert(:,1)';
    candidateqn_n=gas_matrix_vert(:,2)';
    qn_prob=zeros(1,birthn);
    for i=1:birthn
    qn_prob(i)=(1-qn_n_new(i));
    end
    
    logproposal=log((1/prod(qn_prob,'all'))^(Nt+birthn-2));
    Jacobian=(1-sum(qn_n_new))^Nt;
    logJacobian=log(Jacobian);
    %candidatet0=abs(candidatet0);
    %idxt0=find(t0_new>tau);
    %candidatet0(idxt0)=2*tau-candidatet0(idxt0);
    
    

    
else % death step
    Nt=length(currentt0);
    deathn=randperm(2,1);
    
    deathidx=randperm(Nt,deathn);
    candidatet0=currentt0;
    candidatet0(deathidx)=[];
    candidateqn_n=currentqn_n;
    candidateqn_n(deathidx)=[];    
    qn_n_delete=currentqn_n(deathidx);
    candidateqn_n=candidateqn_n/sum(candidateqn_n);
    qn_prob=zeros(1,deathn);
    for i=1:deathn
    qn_prob(i)=(1-qn_n_delete(i));
    end
    logproposal=log((prod(qn_prob,'all'))^(Nt-deathn-2));
    Jacobian=(1/(1-sum(qn_n_delete)))^(Nt-deathn);
    logJacobian=log(Jacobian);
end    


end