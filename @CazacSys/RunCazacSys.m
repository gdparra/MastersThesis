function [varargout] = RunCazacSys(varargin)
dbstop if error

cazacj=varargin{1};
indexval=varargin{2};
TowInd=varargin{3};
DelayToTarget=varargin{4};
CSymbols_PT=varargin{5};
PowerGain=varargin{6};
bind1=varargin{7};
Intbaseband=varargin{8};
h00N=varargin{9};
cazacj.symbind= indexval;
cazacj=freq2time(cazacj,TowInd,DelayToTarget,CSymbols_PT,PowerGain,bind1,Intbaseband,h00N);
varargout{1}=cazacj;




