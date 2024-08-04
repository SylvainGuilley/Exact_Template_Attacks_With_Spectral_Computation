# Context

The so-called template attacks (TA) is one of the optimal Side-Channel Analysis (SCA). In the scientific literature, several optimisations of its implementation are introduced using the so called coalescence, thanks to the Law of Large Numbers (LLN), then using a spectral computation.
Their inconvenient, is that the LLN is just an asymptotic approximation. So it not leads to an exact Template Attack, especially for a small number of traces.

In this paper, we introduce a way of calculating the TA exactly and with the same computational complexity (using the spectral approach), without using the LLN, regardless of the number of messages.
For the experimental validation of this approach, we used the ANSSI SCA Database (ASCAD), with different number of messages and different amount of samples per trace. Recall that this dataset concerns a software implementation of AES-128 bits, running on an ATMEGA-8515 microprocessor.

# Content

In this repository, please find:

- `main_Template_coalescence.c`: previous art approach
- `main.c`: new approach
