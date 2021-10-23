# Note

```bash
conda install -c conda-forge pythonocc-core
pip install PyMieScatt
pip install numdifftools
```

<https://pymiescatt.readthedocs.io/en/latest/examples.html>

## Reflection and Refraction

Snell Raw
$$ n_1\sin(\theta_1) = n_2\sin(\theta_2) $$

- p-pol(TM), s-pol(TE)
- Reflection Index $r_{12}$ define
- Transmit Index $t_{12}$ define

$$ r_{12}^s = \frac{n_1\cos(\theta_1) - n_2\cos(\theta_2)}{n_1\cos(\theta_1) + n_2\cos(\theta_2)} = \frac{k_{1z} - k_{2z}}{k_{1z} + k_{2z}}$$
$$ t_{12}^s = \frac{2n_1\cos(\theta_1)}{n_1\cos(\theta_1) + n_2\cos(\theta_2)} = \frac{2k_{1z}}{k_{1z} + k_{2z}}$$
$$ r_{12}^p = \frac{n_2\cos(\theta_1) - n_1\cos(\theta_2)}{n_2\cos(\theta_1) + n_1\cos(\theta_2)} = \frac{\epsilon_2k_{1z} - \epsilon_1k_{2z}}{\epsilon_2k_{1z} + \epsilon_1k_{2z}}$$
$$ t_{12}^p = \frac{2n_1\cos(\theta_1)}{n_1\cos(\theta_1) + n_2\cos(\theta_2)} = \frac{2n_1n_2k_{1z}}{\epsilon_2k_{1z} + \epsilon_1k_{2z}}$$

Brewster Angle
$$ \tan(\theta_{\Beta}) = \frac{n_2}{n_1} $$

Brewster Angleでは反射率が0、透過率1となる。
高出力レーザーの透過窓は反射光によりレーザー内部が損傷しないように窓がBrewster Angleになるように作られる。

s-polとp-polでは、s-polが常に煩瑣率が高い。
臨界角(Brewster Angle)以降は反射率が常に1になる。
媒質2(入射される側の材質)では消衰波(Evanescent Wave)が存在する。

侵入長
$$ z_d = \frac{\lambda}{2\pi\sqrt{n_1^2\sin^2(\theta_1) - n_2^2}} $$

薄膜の多重反射と透過

3層のは媒質(2層目の媒質の厚みを$d_2$とする)

$$ r_{13} = \frac{r_{12} + r_{23} \exp(2k_{2z}d_2i)}{1 + r_{23}r_{12}\exp(2k_{2z}d_2i)} $$
$$ t_{13} = \frac{t_{12}t_{23}\exp(2k_{2z}d_2i)}{1 + r_{23}r_{12}\exp(2k_{2z}d_2i)} $$
$$ k_{2z} = \frac{2\pi}{\lambda} n_2 \cos(\theta_2) $$

## Propagation matrix of isotropic medium

$$ G = M_{f(f-1)} \Phi_{f-1} ... M_{32} \Phi_2 M_{21} $$

$$ M_{ij}^{P} = \frac{1}{2n_{i}n_{j}k_{iz}} \begin{vmatrix}
    n_i^2k_{iz} + n_j^2k_{iz} & n_i^2k_{jz} - n_j^2k_{iz} \\
    n_i^2k_{iz} - n_j^2k_{iz} & n_i^2k_{jz} + n_j^2k_{iz}
\end{vmatrix} $$

$$ M_{ij}^{S} = \frac{1}{2n_{i}n_{j}k_{iz}} \begin{vmatrix}
    k_{iz} + k_{jz} & k_{iz} - k_{jz} \\
    k_{iz} - k_{jz} & k_{iz} + k_{jz}
\end{vmatrix} $$

$$ \Phi_{i} = \begin{vmatrix}
    \exp(ik_{iz}d_i) & 0 \\
    0 & \exp(-ik_{iz}d_i)
\end{vmatrix} $$

$$ M_{21}^{P} = \frac{1}{2n_{2}n_{1}k_{2z}} \begin{vmatrix}
    n_2^2k_{2z} + n_1^2k_{2z} & n_2^2k_{1z} - n_1^2k_{2z} \\
    n_2^2k_{2z} - n_1^2k_{2z} & n_2^2k_{1z} + n_1^2k_{2z}
\end{vmatrix} $$

$$ M_{21}^{S} = \frac{1}{2n_{2}n_{1}k_{1z}} \begin{vmatrix}
    k_{2z} + k_{1z} & k_{2z} - k_{1z} \\
    k_{2z} - k_{1z} & k_{2z} + k_{1z}
\end{vmatrix} $$
