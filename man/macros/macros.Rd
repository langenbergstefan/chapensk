% Macros for package chapensk

% subscript
\newcommand{\sub}{\ifelse{latex}{\out{$_{#1}$}}{\ifelse{html}{\out{<sub>#1</sub>}}{#1}}}

% greek letters
\newcommand{\Rzeta}{\ifelse{latex}{\out{$\zeta$}}{\ifelse{html}{\out{&zeta;}}{\out{ζ}}}}
\newcommand{\Rsigma}{\ifelse{latex}{\out{$\sigma$}}{\ifelse{html}{\out{&sigma;}}{\out{σ}}}}
\newcommand{\Rvarepsilon}{\ifelse{latex}{\out{$\varepsilon$}}{\ifelse{html}{\out{&epsilon;}}{\out{ε}}}}
\newcommand{\Reta}{\ifelse{latex}{\out{$\eta$}}{\ifelse{html}{\out{&eta;}}{\out{η}}}}
\newcommand{\RTheta}{\ifelse{latex}{\out{$\Theta$}}{\ifelse{html}{\out{&Theta;}}{\out{Θ}}}}
\newcommand{\ROmega}{\ifelse{latex}{\out{$\Omega$}}{\ifelse{html}{\out{&Omega;}}{\out{Ω}}}}

% emph
\newcommand{\Remph}{\ifelse{latex,html}{\emph{#1}}{#1}}