# Note

```bash
conda install -c conda pythonocc-core
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
