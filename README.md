# MTT ICM

Humble C++ lib to compute ICM:
- exactly
- with a Monte-Carlo evaluation for larger numbers

with a big focus on ***performance***.


Extra features:
- evaluate the statistical guarantee along your MC computation
- Python bindings using Boost Python

## History

This code is part of a larger freestyle experiment around valuation functions in poker, that is described in these (french) posts:
- [Poker : MTT et ICM #1 - La Question](https://pittscraft.com/posts/poker_mtt_icm_question/)
- [Poker : MTT et ICM #2 - Le Calcul Brutal](https://pittscraft.com/posts/poker_mtt_icm_calcul/)
- [Poker : MTT et ICM #3 - Méthode de Monte-Carlo](https://pittscraft.com/posts/poker_mtt_icm_monte_carlo/)
- [Poker : MTT et ICM #4 - Perf Tuning avec du Deep Learning](https://pittscraft.com/posts/poker_mtt_icm_deep_learning/)

(hopefully more to come)

As experimented in post #4, training a simple NN on MC-generated data may provide better performance (hundreds of times faster).
The related Python code will come in a separate repo a priori.

## Future

On my own I'll improve and extend this codebase only on need. If you want to contribute, feel free to submit any pull request or directly contact me via email: [pierre@pittscraft.com](mailto:pierre@pittscraft.com).

## Project layout

```
|
|- include    # Public headers
|- src        # Private headers and implementation
|- test       # Boost tests
```

Explore the include folder and you'll find the properly documented few main functions.

## Build

Requires:
- [Boost 1.75.0](https://www.boost.org/users/history/version_1_75_0.html) 
- [Boost Python](https://www.boost.org/doc/libs/1_75_0/libs/python/doc/html/building/installing_boost_python_on_your_.html)
- [SFMT](https://github.com/heavywatal/sfmt-class)

Some additional information is available in [CMakeLists.txt](./CMakeLists.txt).

## License

Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License.

You can use the code freely for any non-commercial use as long as you propagate its license. See [the license file](./license.txt).
