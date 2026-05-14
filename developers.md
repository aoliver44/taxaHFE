![test workflow](https://github.com/aoliver44/taxaHFE/actions/workflows/test.yml/badge.svg)
# **Development** <a><img src='pictures/logo.png' align="right" height="150" /></a>
<a><img src='pictures/dietML_logo_Nov18.png' align="right" height="120" /></a>

The TaxaHFE repo provides an RStudio server image to aid in devlopment. To launch the server use the following command:

```
docker compose up -d
```

To access the server navigate to http://localhost:8787, the username and password will be `rstudio`.

*Note: This will use the remotely build `aoliver44/taxa_hfe_rstudio` by default, if a local build is required, add the `--build` flag to the `docker compose` command above.*

### Restarting/stoppinp

The server can be restarted using:
```
docker compose restart
```

The server can be stopped using:
```
docker compose down
```