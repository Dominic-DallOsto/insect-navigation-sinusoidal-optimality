# Sinusoidal Optimality Paper

Idea: Show that the a sinusoidal activity profile in a circulant ring network has some nice properties - noise resistant and self-stabilising with a Hebbian learning rule.

## Noise Resilience

### Non-negative activity

* [Noise resilience of using fewer Fourier modes](sinusoidal%20noise%20rejection.ipynb)
* [Trying the same with a spiking network](sinusoidal%20noise%20rejection%20spiking.ipynb)

### Zero mean activity

* [Noise resilience of using fewer Fourier modes](sinusoidal%20noise%20rejection%20negative%20activity.ipynb) - Figure 1

## Oja's rule

* [Oja's rule causes the weights to converge to a stable sinusoidal profile](sinusoidal%20Oja's%20rule%20negative%20activity.ipynb)

## Hebbian learning

* [Hebbian learning causes the weights to converge to a stable sinusoidal profile](sinusoidal%20Hebbian%20learning%20weights.ipynb)
* [Phase portraits for the activity and the weights](phase%20portaits.ipynb)
* [Updating both the activity and weights simultaneously](sinusoidal%20updating%20activity%20and%20weights.ipynb) - Figure 5

### Zero mean activity

* [Hebbian learning causes the weights to converge to a stable sinusoidal profile](sinusoidal%20Hebbian%20learning%20weights%20negative%20activity.ipynb)

## Comparison to experimental data

* [comparing the 4 population networks from Ioannis's paper to the abstract network we predict](insect%20network%20graph%20analysis.ipynb) - Figures 3 and 4
* [deriving the 4 population connectivity network from Janelia's connectome data](fruit_fly_data_analysis.ipynb)

## Multiple harmonic encodings

* [What does it look like to encode with multiple harmonics at once](encoding%20with%20multiple%20harmonics.ipynb) - Figure 2
* [How is an encoding with multiple harmonics affected by noise](noise%20affected%20multiple%20harmonics%20encoding.ipynb) - Figures 6 and 7


## Plot styling

* use seaborn-notebook style
* full linewidth - make the plot 10 inches wide so the text is nicely readable in the paper
* use TeX rendering
* nice neuron activity plots - `'o--', markersize=6, linewidth=1`

For poster / when you need a different font in the figure but also LaTeX rendering - save the plot as .pgf and it'll render with the right font in LaTeX.
