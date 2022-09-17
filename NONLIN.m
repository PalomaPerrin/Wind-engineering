function f=NONLIN(phiy,phiz,phit,qn0,Cd,Cl,Cm,rho,B,L,alpha,K_m,U)      
        
            tetha=phit*qn0;% Mi aspetto una 243 per 1
        
                for i=1:243
                    Fy(i)=0.5*rho*(U^2)*(B)*interp1(alpha,Cd,tetha(i)); % U è un valore solo
                    Fz(i)=0.5*rho*(U^2)*(B)*interp1(alpha,Cl,tetha(i));
                    M(i)=0.5*rho*(U^2)*(B^2)*L*interp1(alpha,Cm,tetha(i));
                end
        
            f=K_m*qn0-(phiy'*Fy'+phiz'*Fz'+phit'*M');       
         
        
  
  
  
  
      