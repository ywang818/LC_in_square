function M=find_prc_INI(model)


zinit_matrix=eye(2);
M=[];

for i=1:2
    model.find_prc(zinit_matrix(i,:));
    M=[M; model.prc(end,:)]; %Monodromy matrix
end

end