function M=find_prc_monodromy(model)
% Utility function for Fig5B_prc_plot.m
% This function finds the monodromy matrix associated with the iPRC

zinit_matrix=eye(2);
M=[];

for i=1:2
    model.find_prc(zinit_matrix(i,:));
    M=[M; model.prc(end,:)]; %Monodromy matrix
end

end