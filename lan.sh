brew services restart shopify/shopify/toxiproxy
toxiproxy-cli create SLOW --listen localhost:55555 -u localhost:55556
toxiproxy-cli toxic add SLOW -t latency -a latency=2
toxiproxy-cli toxic add SLOW -t bandwidth -a rate=125000
