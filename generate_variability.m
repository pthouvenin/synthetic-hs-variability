function m_var = generate_variability(m,L,coeffvar)

L_break = floor(L/2 + floor(rand*L/3));
Xi = unifrnd(1-coeffvar/2,1+coeffvar/2,3,1);

y1 = Xi(1);
y2 = Xi(2);
y3 = Xi(3);

vect_mult = [linspace(y2,y1,L_break), linspace(y1,y3,L-L_break)]'; 
m_var = vect_mult.*m;

end