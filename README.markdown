# LogStack

LogStack is the state-of-the-art technique for Stacked Garbling.
Stacked Garbling is a Garbled Circuit improvement that accelerates the handling of conditional branches.
For more information, have a look at our paper ([Heath and Kolesnikov, Eurocrypt 2021](https://eprint.iacr.org/2021/531)).

## Building and Running

First, you'll need to build and install the emp toolkit and the emp OT library: https://github.com/emp-toolkit

- `git submodule init`
- `git submodule update`
- `cd p`
- `make`


For benchmarking, we throttle the network using [Toxiproxy](https://github.com/Shopify/toxiproxy).
Make sure to install that the software so you can run `toxiproxy-cli` on the command line.
