%Created by Michael Leznik
%Distributed under GNU GENERAL PUBLIC LICENSE
%When using this code please refer/cite  
%Estimating Invariant Principal Components Using Diagonal Regression.
%Fits orthognal planes using the least volumes criteria
%Performs nonlinear constrained minimization, using MATLAB
%procedure fmincon
%startPoints - starting values have to be passed in array
%dataSet passed as matrix


function volcomps = volumes(dataSet, startPoints)
    global tempAnswers;
    temp=[];
    meanArr=[];
    [row, col]=size(dataSet);
    for n=1:col %centers all the data by subtracting respective mean
        meanArr(n)= mean(dataSet(:,n));
        for r=1:row
            dataSet(r,n)=dataSet(r,n)-meanArr(n);
        end    
    end
    if (nargin >1 && (length(startPoints) == col))
        a0=startPoints; % starting values specified by user
    else
        for v=1:col %starting values given by user, all values set to 1.0
            a0(v)=1.0;
        end    
    end
    for run=1:col 
        [temp, fval]=mvolfunc(dataSet, run, a0);
        volcomps(:,run)=temp; %values to the function return
        tempAnswers(:,run)=temp;%assigns answers to global array in order  
                                                     % to use them for constraints
    end
    clear global answers; % clear global array after function finishes run
end
function [a, fval]=mvolfunc(dataSet, run, a0)
    options = optimset('LargeScale','off','MaxIter', 10000,...
                        'TolCon',1.e-7,'TolX',1.e-7,'Display','iter',...
                        'MaxFunEvals', 10000,'TolFun',1.e-7);
    [row, col]=size(dataSet); 
    [a, fval] = fmincon(@volumefunc,a0,[],[],[],[],[],[],...
                        @constrfun,options); % calls constrained 
                                              %minimization function
    
    function f = volumefunc(a)
        % Nonlinear objective function
        Rfunc = 0;
        temp = [];
        temp1=[];
        temp1=1;      
        for l = 1:row
            temp=0;
            for k=1:col
                 temp= temp + (dataSet(l,k)*a(k));
            end
            Rfunc= Rfunc + temp^(2*col);
        end
        for v=1:col  %create denominator of objective function
            temp1 = temp1 * a(v); 
        end    
        f=Rfunc/temp1^row;
        f=Rfunc;
    end    
    function [c,ceq] = constrfun(a)% function creates constraints 
        global tempAnswers;
        temp1=[];
        temp2=[];
        temp3=[];
        temp4=[];
        temp4=1;
        c=[];  %Nonlinear inequality constraints
        ceq=[];%Nonlinear equality constraints
        for v=1:col
            temp4 = temp4 + a(v)^2; %sum of squares equals 1 constraint
            %temp4 = temp4 * a(v); % product equals 1 constraint
        end    
        ceq(1)=temp4-1;
        for s=run:-1:2 % runs backwards and creates orthogonality 
                       % constraints. To make sure that each plane 
                       % is orthogonal with the previous planes
            temp1=0;
            temp2=0;
            temp3=0;
            for l = 1:row
                for k=1:col
                    temp1 = temp1 + (dataSet(l,k)*a(k));
                    temp2 = temp2 + (dataSet(l,k)*tempAnswers(k,s-1));
                end
                temp3 = temp3 + (temp1*temp2); % sums products
                temp1=0;
                temp2=0;
            end
           ceq(s)=temp3;
        end    
    end  
end