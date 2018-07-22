function A = smooth_abundance(sigma2,H,W,P,cutoff,position)
%-Génération de matrices d'abondance synthétiques comportant une forte 
% corrélation spatiale (smoothness).
%%-----------------------------------------------------------------------%%
%-Arguments :
%     H    : hauteur de la carte d'abondance générée;
%     W    : largeur   --------------------------   ;
%     P    : nombre d'endmembers considérés;
%   cutoff : somme max. des abondances par pixel (absence de pixel pur);
%   sigma2 : règle la rapidité à laquelle les abondances décroissent;
%   method : positions des pixels purs aléatoires ou fixées
%
%-Sortie : 
%   A : matrice de format (P|H*W) contenant les cartes d'abondances
%       générées pour chacun des endmembers.
%%-----------------------------------------------------------------------%%
A = zeros(H,W,P);
for k = 1:P-1
    for i = 1:H
        for j = 1:W
            A(i,j,k) = cutoff*exp( -(norm([i,j] - position(k,:),2)^2)/(2*sigma2));
        end
    end
end

A(:,:,P) = ones(H,W) - sum(A(:,:,1:P-1),3);
id  = (A > cutoff);
A(id) = cutoff;
id = (A <= 0);
A(id) = 0;

s = sum(A,3);
A = A./s(:,:,ones(P,1));

A = (reshape(permute(A,[2 1 3]),H*W,P))'; % remise des éléments dans l'ordre (lié au choix de l'ordre lexicographique)

end