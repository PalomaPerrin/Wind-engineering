function PHI=PHI_gen(fisez,modo)
    
        fisez=permute(fisez,[3,1,2]);
        
             for i=1:243
                 for j=1:3
                     PHI(i,j)=fisez(i,j,modo);
                 end
             end