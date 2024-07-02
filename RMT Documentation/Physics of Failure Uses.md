# Physics of Failure Model and Use Table

| Category | Model | Equation | Parameters | Physics of Failure Applications | RMT Tool(s) |
| -------- | ----- | -------- | ---------- | ------------------------------- | -------- |
| **Stress-Life** | S-N | $$S_{a} = A N_{f}^b \text{ or } N_{f} = \left(\frac{S_{a}}{A}\right)^{\frac{1}{b}}$$ | $$A, b$$ | Low cycle fatigue |      |
| **Strain-Life** | Coffin-Manson Relationship | $$\epsilon_a = \frac{\sigma_f'}{E} (2 N_f)^b + \epsilon_f' (2 N_f)^c$$ | $$\epsilon_a$$ Total strain amplitude, $$E$$ Modulus of elasticity, $$\sigma_f'$$ | High cycle fatigue |      |
|              | Morrow Mean-Stress Correction | $$\epsilon_a = \frac{\sigma_f' - \sigma_m}{E} (2 N_f)^b + \epsilon_f' (2 N_f)^c$$ |    |      |      |
|              | Modified Morrow Mean-Stress Correction | $$\epsilon_a = \frac{\sigma_f' - \sigma_m}{E} (2 N_f)^b + \epsilon_f' (\frac{\sigma_f' - \sigma_m}{\sigma_f'})^{\frac{c}{b}} (2 N_f)^c$$ |    |      |      |
|              | Smith-Watson-Topper Mean-Stress Correction | $$\sigma_{max} \epsilon_a = (\sigma_a + \sigma_m) \epsilon_a = \frac{\sigma_f'^{2}}{E} (2 N_f)^{2b} + \sigma_f'\epsilon_f' (2 N_f)^{b+c}$$ |    |      |      |
|              | Walker Mean-Stress Correction | $$\epsilon_a = \left(\frac{\sigma_f'}{E}\right) \left(\frac{1 - R}{2}\right)^{1-\gamma} (2 N_f)^b + \epsilon_f' \left(\frac{1 - R}{2}\right)^{\frac{c}{b}\left(1-\gamma\right)} (2 N_f)^c$$ |    |      |      |
| **Variable Amplitude Loading** | Palmgren's Miner's Linear Damage Rule | $$B_f\left(\sum_{i=1}\frac{n_i}{N_{f,i}} \right)_{\text{one repetition}} = 1$$ |  | Low/High cycle fatigue |      |
|  |  | $${D_i = B_f\left(\frac{n_i}{N_{f,i}}\right)}$$ |  |  |      |
|  | Kwofie-Rahbar Non-linear Damage Rule | $$B_f\left[\sum_{i=1}\frac{n_i}{N_{f,i}} \frac{\ln\left(N_{f,i}\right)}{\ln\left(N_{f,1}\right)} \right]_{\text{one repetition}} = 1$$ |  |  |      |
|  |  | $$D_i = B_f\left[\frac{n_i}{N_{f,i}} \frac{\ln\left(N_{f,i}\right)}{\ln\left(N_{f,1}\right)}\right]$$ |  |  |      |
| **Fracture Mechanics** | Threshold Model | $$\frac{da}{dN} = A \left(\Delta K - \Delta K_{th}\right)^p$$ |  | Region I |      |
|  | Paris-Edrogan Model | $$\frac{da}{dN} = C \left(\Delta K\right)^m$$ | $$C, m$$ | Region II |      |
|  | Walker Model | $$\frac{da}{dN} = C \left[\frac{\Delta K}{(1-R)^{1-\gamma}}\right]^m$$ | $$C, m$$ | Region II |      |
|  | Foreman Model | $$\frac{da}{dN} = \frac{C_{III} \left(\Delta K\right)^{m_{III}}}{\left(1-R\right) K_{ic}-\Delta K}$$ | $$C, m$$ | Region III |      |
|  | Mechanistic or Exponential Model | $$\frac{da}{dN} = ma \text{ or } a(N) = a_0 \exp{mN}$$ | $$C, m$$ | Region I, II, & III |      |
|  | McEvily-Groegrer Model | $$\frac{da}{dN} = A \left(\Delta K - \Delta K_{th}\right)^2 \left[1+\frac{\Delta K}{K_{ic}-K_{max}}\right]$$ | $$C, m$$ | Region I, II, & III |      |
|  | NASGRO Model | $$\frac{da}{dN} = C \left[\left(\frac{1-f}{1-R}\right)\Delta K\right]^n \frac{\left(1-\frac{\Delta K_{th}}{\Delta K}\right)^p}{\left(1-\frac{K_{max}}{K_{ic}}\right)^q}$$ | $$C, m$$ | Region I, II, & III |      |
| **Wear** | Archard's Wear Law | $$W=k\frac{NL}{H}$$ | $$N, L, H$$ |  |      |
