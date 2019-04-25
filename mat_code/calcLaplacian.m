function [eigVector,eigValue] = calcLaplacian(L)%,cmd
% figure;
 % Calculating Laplacian matrices
%add by shiqun，自动获取顶点数目，不用传入参数
[eigVector, eigValue]  = eig(L); %describe by shiqun,eigs是matlab自带的做特征分解的函数?
eigValue = diag(eigValue) - 1;
eigValue = eigValue(end:-1:1);eigValue = max(eigValue(:),0); % the first element due to num. error is negative: -1e-15
eigVector = eigVector( :,end:-1:1 );
% figure;   
% imagesc (eigValue);
% figure;
% imagesc (eigVector);,struct('disp', 0)
end