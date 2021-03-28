function[QFinal,NFinal,exitflagfzero,exitflagRelax] = NewRaphToFindQtop3(Nbottom,S,psi,ds,s,eps,r,delta,omega,NL,QLast,DPing,PlotRelax)
% Qtop is the guess at the top Q
% Nbottom is the bottom BC on NR
% S is the tunnel cross-sectional area profile
% NC is the cavity effective pressure
% psi is the basic hydraulic gradient
% ds is the space step
% s is the space vector
% omega is the input to the tunnnel per unit length per unit time
% kapR is the connectivity between the two drainage systems
% NL is the lake effective pressure
% QLast os the Q from thje last time step to use in the initial guess for this one
%
% See section 2.25 of Kinglake (2013) for details of the method:
% http://etheses.whiterose.ac.uk/4630/1/J_Kingslake_ThesisFINAL.pdf
%
% Written by J. Kingslake, 2011, Department of Geography, University of
% Sheffield, UK.




% For the inital guess that you need to start the N-R method use the
% previous timesteps Q at the top
% [QtopFinal,TotIter,exitflag,QFinal,NFinal] = NewtonRaphson_RelaxMeth_2(@NtopFromQtopGuess,QLast(1,1));
%
%     warning 'Problem with the Newton-Raphson Method'
%     if there has been a problem try using the built in root finding function fzero
% if exitflag~=1
[QtopFinal,fval,exitflagfzero,output] = fzero(@NtopFromQtopGuess,QLast(1,1));
[NtopMinusNL,Ntop,exitflagRelax,QFinal,NFinal] = NtopFromQtopGuess(QtopFinal);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%  This version doesnt deal with coupling between the channel and a
%%%%%%%  cavity system!!!!!!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function[NtopMinusNL,Ntop,exitflag,QFinal,NFinal] = NtopFromQtopGuess(QtopGuess);
        %%%%% RELAXATION METHOD %%%%%%%%%
        % See section 2.25 of Kinglake (2013)
        % http://etheses.whiterose.ac.uk/4630/1/J_Kingslake_ThesisFINAL.pdf
        % the output Ntop is the effective pressure at the top of the tunnel which
        % results from the given Qtop, S, psi etc.
        % exitflag is a variable which indicates how the function finished
        % =1 if the solution relaxed to a 'steady-state' solution
        % =-1 if it reached the max number of iterations without converging
        % =0 neither happened, exitflag wasnt defined after the initiatioj of this function
        
        % isnt the same time as in the main function. It is just a way of allowing
        % the solution to relax to the correct solution.
        
        T_RM = 50000;
        %dt_RM = 0.001;        % use this one for Gorner
        dt_RM = 0.005;          % choose this: lower increase instability of solution noise, whereas higher makes this method converge faster.
        % t = 0:dt_RM:T_RM;
        % DPing = 0.05;
        exitflag = 0;
        Ls = length(s);
        CS =zeros(1,Ls-1);
        % tolerance, when the top N changes less than this in one timestep the
        % function has found a 'steady-state' solution
        Tol_RM = 1e-9;         % choose this: lower reduces numerical noise, whereas higher makes this method converge faster.
        
        Q= ones(2,Ls);
        N = zeros(2,Ls) + Nbottom;
        ddtN = zeros(1,Ls);
        
        % pre-do the S^(8/3)
        S_83 = S.^(8/3);
        % and dt_RM/delta
        A = ds/delta;
        
        % pre-allocated constants to speed things up
        C1 = A./S_83(1,2:Ls);
        C2 = A.*psi(1,2:Ls);
        C3 = eps*(r-1)./S_83(1,2:Ls);
        C4 = eps*S(1,2:Ls);
        
        
        
        % inital guess at Q has the current guess at Q at the top and the last
        % timesteps Q profile for the rest of the tunnel
        Q(1,:) = [QtopGuess QLast(2:Ls)];
        % Q(1,:) = QtopGuess*ones(1,Ls);
        err = 0;
        % NCShort = NC(1,2:Ls);
        for i = 1:T_RM,
            %      T_RM-i
            QShort = Q(1,2:Ls);
            
            %     Q_sqr = realpow(QShort,2);
            Q_sqr = QShort.*QShort;
            vec1 = -(C1.*Q_sqr.*sign(QShort)- C2);
            
            % the two lines below are a replacement for the loop which runs over
            % the space domain. It is quicker to do it in this new way!!
            CS = cumsum(vec1(Ls-1:-1:1));
            N(2,Ls-1:-1:1) = Nbottom + CS;
            
            %     for j = Ls:-1:2,
            %         N(2,j-1) = N(2,j) + vec1(j-1); % the use of (j-1) in the vec1 index is because vec1 is already only Ls-1 elements long and the jth element of vec1 actually refers to the j+1th position in space
            %     end
            NShort =  N(2,2:Ls);
            %     NShort_sqr = realpow(NShort,2);
            SourceTerms =  omega;
            ddsQ = (QShort - Q(1,1:Ls-1))/ds;
            
            v1= -ddsQ;
            v2= C3.*abs(QShort).*Q_sqr;
            %     v3= C4.* realpow(NShort,3);
            %     v3= C4.* NShort_sqr .* NShort;
            v3= C4.* NShort .* NShort.* NShort;
            %      v3= C4.* NShort.^3;
            v4= SourceTerms;
            % 'time'-derivative of N
            
            if DPing ==0
                vec2 = v1+v2+v3+v4;
            else
                if i>1
                    ddtN = (N(2,:)-N(1,:))/dt_RM;
                end
                vec2 = v1+v2+v3+v4+DPing.*S(1,2:Ls).*ddtN(1,2:Ls);
            end
            %     vec2 = -ddsQ + C3.*abs(QShort).*Q_sqr + C4.* N(2,2:Ls).^3 + SourceTerms ;
            
            Q(2,2:Ls) = QShort + dt_RM * vec2;
            Q(2,1) = Q(1,1);
            %  Q(2,:) = [Q(1,1) QShort + dt_RM * vec2];
            
            N(2,Ls) = N(1,Ls);
            
            if PlotRelax == 1 && rem(i,100)==0
                figure(1)
                %    hold on
                plot(s,Q(2,:)); title('Q')
                figure(2);
                plot(s,N(2,:)); title('N')
                %     hold off
            end
            if i > 2 && rem(i,10)==0
                %         err = abs(N(2,:) - N(1,:))/N(2,:);
                err1 = max(abs(N(2,:) - N(1,:))./N(2,:));
                err2 = max(vec2);
                % maxerr = max(abs(ddtN));
                %         err = max(abs(Q(2,:)-Q(1,:))) + max(abs(N(2,:)-N(1,:)));
                if abs(err1+err2) < Tol_RM
                    IterToConverge = i;
                    exitflag =1;
                    Ntop = N(2,1);
                    NtopMinusNL = N(2,1) - NL;
                    NFinal = N(2,:);
                    QFinal = Q(2,:);
                    
                    return
                end
            end
            % loop values round
            Q(1,:) = Q(2,:);
            N(1,:) = N(2,:);
        end
        warning 't reached T_RM without converging!! and NaN returned for the value of Ntop'
        save 'InnerSave of fail'
        exitflag = -1;
        Ntop = NaN;
        NtopMinusNL = NaN;
        Nfinal = NaN;
        QFinal = NaN;
        
    end
end