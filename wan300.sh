brew services restart shopify/shopify/toxiproxy
toxiproxy-cli create SLOW --listen localhost:55555 -u localhost:55556
# toxiproxy-cli toxic add SLOW -t latency -a latency=100
toxiproxy-cli toxic add SLOW -t bandwidth -a rate=37500
# toxiproxy-cli toxic add SLOW -t bandwidth -a rate=6250