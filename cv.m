function VV = cv(yita,V1,V2,V3,V4,V)
VV=(yita-1)*V1+yita*V2/2+yita*(1-yita)*V3/6+yita*(1-yita)*(2-yita)*V4/24+V;
end
