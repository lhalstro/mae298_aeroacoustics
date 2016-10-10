
%read audio
[y,Fs] = audioread('Data/Boom_F1B2_6.wav');
disp(max(y))
dlmwrite('Data/MatlabData.dat', y, ' ')

%play raw sound (non-normalized)
pobj = audioplayer(y, Fs);
play(pobj);
pause(1);

%audioplayer in matlab
fs = 44100;
%normalize sound so you can hear peaks better
pvolts = y / max(abs(y));
pvolts = [pvolts];
ppobj = audioplayer(pvolts, fs);
play(ppobj);