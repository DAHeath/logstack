#include <emp-ot/iknp.h>
#include "protocol.h"


void send(const Label& msg, emp::NetIO& channel) {
  channel.send_block((const emp::block*)&msg, 1);
}


void send(std::span<const Label> msg, emp::NetIO& channel) {
  channel.send_block((const emp::block*)&msg[0], msg.size());
}


Label recv(emp::NetIO& channel) {
  emp::block buf;
  channel.recv_block((emp::block*)&buf, 1);
  return *(Label*)&buf;
}


std::vector<Label> recv(emp::NetIO& channel, std::size_t n) {
  std::vector<Label> out(n);
  channel.recv_block((emp::block*)out.data(), n);
  return out;
}


void semihonestGen(const Circuit& c, emp::NetIO& channel) {
  PRG prg;

  const auto fixedKey = prg();

  const PRF f { fixedKey };

  const auto [m, i] = garble(f, c, prg);
  const auto [e, d] = i;

  send(fixedKey, channel);

  send(m, channel);


  const auto n = c.nInp;
  std::vector<emp::block> data0(n);
  std::vector<emp::block> data1(n);
  for (std::size_t i = 0; i < n; ++i) {
    const auto l = e.zeros[i];
    const auto r = e.zeros[i] ^ e.delta;

    data0[i] = *(emp::block*)&l;
    data1[i] = *(emp::block*)&r;
  }


  emp::IKNP ot(&channel);
  ot.send(data0.data(), data1.data(), n);

  const auto Y = recv(channel, c.nOut);

  /* for (std::size_t i = 0; i < c.nOut; ++i) { */
  /*   if (Y[i] == d.zeros[i]) { */
  /*     std::cout << '0'; */
  /*   } else if (Y[i] == (d.zeros[i] ^ d.delta)) { */
  /*     std::cout << '1'; */
  /*   } else { */
  /*     /1* std::cerr << "ERROR!\n"; *1/ */
  /*     /1* std::cerr << out[i] << '\n'; *1/ */
  /*     /1* std::cerr << g.outputEncoding.zeros[i] << '\n'; *1/ */
  /*     /1* std::cerr << (g.outputEncoding.zeros[i] ^ delta2) << '\n'; *1/ */
  /*     /1* std::exit(1); *1/ */
  /*   } */
  /* } */
  /* std::cout << '\n'; */
}


void semihonestEval(const Circuit& c, emp::NetIO& channel) {
  const auto fixedKey = recv(channel);
  const PRF f { fixedKey };

  auto m = recv(channel, c.nRow);


  const auto n = c.nInp;
  bool* x = new bool[n];
  std::vector<Label> X(n);
  for (std::size_t i = 0; i < n; ++i) {
    x[i] = 0;
  }

  emp::IKNP ot(&channel);
  ot.recv((emp::block*)X.data(), x, n);

  delete[] x;

  const auto Y = ev(f, c, X, m);

  send(Y, channel);
  channel.flush();
}



void semihonest(const Circuit& c, emp::NetIO& genChannel, emp::NetIO& evalChannel) {
  std::thread th { [&] {
    semihonestGen(c, genChannel);
  }};

  semihonestEval(c, evalChannel);
  th.join();

}
