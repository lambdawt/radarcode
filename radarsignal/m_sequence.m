clear;
close all;
% r=3; gen=[1 0 1 1];
% r=4; gen=[1 0 0 1 1];
% m序列码生成
r=4; gen=[1 0 0 1 1];
N = 2^r - 1;
reg = [zeros(1, r - 1) 1];  % 初始化register(0...1)
% 参数初始化
seq = zeros(1, N);
newReg = zeros(1, r);

seq(1) = reg(r);
for i = 2 : N
    media = gen .* [0, reg];
    for j = 1:(r - 1)
        newReg(j + 1) = reg(j);
    end
    newReg(1) = mod(sum(media), 2);
    seq(i) = newReg(r);
    reg = newReg;
end
seq    
