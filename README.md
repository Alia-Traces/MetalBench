# MetalBench

A GPU compute benchmarking tool for macOS / Metal

## About

MetalBench contains several 3D scenes and utilises GPU compute shaders to render them using path tracing. Some scenes are rendered realtime, others are offline with incremental display update.

The app uses accurate GPU draw timing to calculate a Rays per Second value.

## How to use

Select a GPU and scene, and observe the Average Rays Per Second value at the bottom right. It's best to wait several minutes, as some scenes contain several views which affect rendering speed, and some GPUs become thermally constrained after a short time and will slow down.

