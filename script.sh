#!/bin/bash

time ./casanova & ./get_client & ./put_client

rm server_get_socket server_put_socket