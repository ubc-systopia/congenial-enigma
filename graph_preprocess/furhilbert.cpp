//
// Created by atrostan on 16/09/22.
//
unsigned long long HILLOOP_nanoprog[9][9][4][2] = { {
		                                                    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 0 x 0
		                                                    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 0 x 1
		                                                    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 0 x 2
		                                                    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 0 x 3
		                                                    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 0 x 4
		                                                    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 0 x 5
		                                                    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 0 x 6
		                                                    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 0 x 7
		                                                    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} } }, // 0 x 8
                                                    { { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 1 x 0
		                                                    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 1 x 1
		                                                    { {0x0ull, 0x2ull}, {0x0ull, 0x2ull}, {0x2ull, 0x2ull}, {0x2ull, 0x2ull} }, // 1 x 2
		                                                    { {0x0ull, 0x4ull}, {0x0ull, 0x4ull}, {0x6ull, 0x4ull}, {0x6ull, 0x4ull} }, // 1 x 3
		                                                    { {0x0ull, 0x8ull}, {0x0ull, 0x8ull}, {0xeull, 0x8ull}, {0xeull, 0x8ull} }, // 1 x 4
		                                                    { {0x0ull, 0x10ull}, {0x0ull, 0x10ull}, {0x1eull, 0x10ull}, {0x1eull, 0x10ull} }, // 1 x 5
		                                                    { {0x0ull, 0x20ull}, {0x0ull, 0x20ull}, {0x3eull, 0x20ull}, {0x3eull, 0x20ull} }, // 1 x 6
		                                                    { {0x0ull, 0x40ull}, {0x0ull, 0x40ull}, {0x7eull, 0x40ull}, {0x7eull, 0x40ull} }, // 1 x 7
		                                                    { {0x0ull, 0x80ull}, {0x0ull, 0x80ull}, {0xfeull, 0x80ull}, {0xfeull, 0x80ull} } }, // 1 x 8
                                                    { { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 2 x 0
		                                                    { {0x0ull, 0x3ull}, {0x0ull, 0x3ull}, {0x2ull, 0x3ull}, {0x2ull, 0x3ull} }, // 2 x 1
		                                                    { {0x8ull, 0xdull}, {0x8ull, 0xaull}, {0x6ull, 0xdull}, {0x6ull, 0xaull} }, // 2 x 2
		                                                    { {0x0ull, 0x1ull}, {0x30ull, 0x24ull}, {0x0ull, 0x1ull}, {0xeull, 0x24ull} }, // 2 x 3
		                                                    { {0x88ull, 0xd5ull}, {0xe0ull, 0x88ull}, {0x76ull, 0xd5ull}, {0x1eull, 0x88ull} }, // 2 x 4
		                                                    { {0x0ull, 0x1ull}, {0x3c0ull, 0x210ull}, {0x0ull, 0x1ull}, {0x3eull, 0x210ull} }, // 2 x 5
		                                                    { {0x888ull, 0xd55ull}, {0xf80ull, 0x820ull}, {0x776ull, 0xd55ull}, {0x7eull, 0x820ull} }, // 2 x 6
		                                                    { {0x0ull, 0x1ull}, {0x3f00ull, 0x2040ull}, {0x0ull, 0x1ull}, {0xfeull, 0x2040ull} }, // 2 x 7
		                                                    { {0x8888ull, 0xd555ull}, {0xfe00ull, 0x8080ull}, {0x7776ull, 0xd555ull}, {0x1feull, 0x8080ull} } }, // 2 x 8
                                                    { { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 3 x 0
		                                                    { {0x0ull, 0x7ull}, {0x0ull, 0x7ull}, {0x6ull, 0x7ull}, {0x6ull, 0x7ull} }, // 3 x 1
		                                                    { {0x30ull, 0x3bull}, {0x0ull, 0x1ull}, {0xeull, 0x3bull}, {0x0ull, 0x1ull} }, // 3 x 2
		                                                    { {0x188ull, 0x1caull}, {0x188ull, 0x135ull}, {0x76ull, 0x1caull}, {0x76ull, 0x135ull} }, // 3 x 3
		                                                    { {0x708ull, 0xa8aull}, {0x0ull, 0x1ull}, {0x8f6ull, 0xa8aull}, {0x0ull, 0x1ull} }, // 3 x 4
		                                                    { {0x6230ull, 0x729bull}, {0x3e20ull, 0x68d4ull}, {0x1dceull, 0x729bull}, {0x41deull, 0x68d4ull} }, // 3 x 5
		                                                    { {0x30c30ull, 0x3b6dbull}, {0x0ull, 0x1ull}, {0xf3ceull, 0x3b6dbull}, {0x0ull, 0x1ull} }, // 3 x 6
		                                                    { {0x188c30ull, 0x1ca6dbull}, {0x1be208ull, 0x128d45ull}, {0x773ceull, 0x1ca6dbull}, {0x41df6ull, 0x128d45ull} }, // 3 x 7
		                                                    { {0xc30c30ull, 0xedb6dbull}, {0x0ull, 0x1ull}, {0x3cf3ceull, 0xedb6dbull}, {0x0ull, 0x1ull} } }, // 3 x 8
                                                    { { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 4 x 0
		                                                    { {0x0ull, 0xfull}, {0x0ull, 0xfull}, {0xeull, 0xfull}, {0xeull, 0xfull} }, // 4 x 1
		                                                    { {0xe0ull, 0xf7ull}, {0x88ull, 0xaaull}, {0x1eull, 0xf7ull}, {0x76ull, 0xaaull} }, // 4 x 2
		                                                    { {0x0ull, 0x1ull}, {0x708ull, 0xd75ull}, {0x0ull, 0x1ull}, {0x8f6ull, 0xd75ull} }, // 4 x 3
		                                                    { {0x7888ull, 0xad5aull}, {0x7888ull, 0xd2a5ull}, {0x8776ull, 0xad5aull}, {0x8776ull, 0xd2a5ull} }, // 4 x 4
		                                                    { {0x0ull, 0x1ull}, {0x7c308ull, 0xd1245ull}, {0x0ull, 0x1ull}, {0x83cf6ull, 0xd1245ull} }, // 4 x 5
		                                                    { {0x1f8888ull, 0x88d55aull}, {0x778888ull, 0xd52a55ull}, {0xe07776ull, 0x88d55aull}, {0x887776ull, 0xd52a55ull} }, // 4 x 6
		                                                    { {0x0ull, 0x1ull}, {0x77c3088ull, 0xd512455ull}, {0x0ull, 0x1ull}, {0x883cf76ull, 0xd512455ull} }, // 4 x 7
		                                                    { {0x78887888ull, 0xad5a2d5aull}, {0x77788888ull, 0xd552a555ull}, {0x87778776ull, 0xad5a2d5aull}, {0x88877776ull, 0xd552a555ull} } }, // 4 x 8
                                                    { { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 5 x 0
		                                                    { {0x0ull, 0x1full}, {0x0ull, 0x1full}, {0x1eull, 0x1full}, {0x1eull, 0x1full} }, // 5 x 1
		                                                    { {0x3c0ull, 0x3efull}, {0x0ull, 0x1ull}, {0x3eull, 0x3efull}, {0x0ull, 0x1ull} }, // 5 x 2
		                                                    { {0x3e20ull, 0x572bull}, {0x6230ull, 0x4d64ull}, {0x41deull, 0x572bull}, {0x1dceull, 0x4d64ull} }, // 5 x 3
		                                                    { {0x7c308ull, 0xaedbaull}, {0x0ull, 0x1ull}, {0x83cf6ull, 0xaedbaull}, {0x0ull, 0x1ull} }, // 5 x 4
		                                                    { {0x7e2308ull, 0x12729baull}, {0x7e2308ull, 0x1d8d645ull}, {0x181dcf6ull, 0x12729baull}, {0x181dcf6ull, 0x1d8d645ull} }, // 5 x 5
		                                                    { {0x7f0c308ull, 0x223b6dbaull}, {0x0ull, 0x1ull}, {0x380f3cf6ull, 0x223b6dbaull}, {0x0ull, 0x1ull} }, // 5 x 6
		                                                    { {0x1f88c30e0ull, 0x49ca6db88ull}, {0x1e7e23088ull, 0x76d8d6455ull}, {0x60773cf1eull, 0x49ca6db88ull}, {0x6181dcf76ull, 0x76d8d6455ull} }, // 5 x 7
		                                                    { {0x7c3087c308ull, 0xaedba2edbaull}, {0x0ull, 0x1ull}, {0x83cf783cf6ull, 0xaedba2edbaull}, {0x0ull, 0x1ull} } }, // 5 x 8
                                                    { { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 6 x 0
		                                                    { {0x0ull, 0x3full}, {0x0ull, 0x3full}, {0x3eull, 0x3full}, {0x3eull, 0x3full} }, // 6 x 1
		                                                    { {0xf80ull, 0xfdfull}, {0x888ull, 0xaaaull}, {0x7eull, 0xfdfull}, {0x776ull, 0xaaaull} }, // 6 x 2
		                                                    { {0x0ull, 0x1ull}, {0x30c30ull, 0x24924ull}, {0x0ull, 0x1ull}, {0xf3ceull, 0x24924ull} }, // 6 x 3
		                                                    { {0x778888ull, 0xaad5aaull}, {0x1f8888ull, 0xf72aa5ull}, {0x887776ull, 0xaad5aaull}, {0xe07776ull, 0xf72aa5ull} }, // 6 x 4
		                                                    { {0x0ull, 0x1ull}, {0x7f0c308ull, 0x3dc49245ull}, {0x0ull, 0x1ull}, {0x380f3cf6ull, 0x3dc49245ull} }, // 6 x 5
		                                                    { {0x1f7888e08ull, 0x88ad5a77aull}, {0x1f7888e08ull, 0xf752a5885ull}, {0xe087771f6ull, 0x88ad5a77aull}, {0xe087771f6ull, 0xf752a5885ull} }, // 6 x 6
		                                                    { {0x0ull, 0x1ull}, {0x21ddf0c3088ull, 0x2b568492455ull}, {0x0ull, 0x1ull}, {0x1e220f3cf76ull, 0x2b568492455ull} }, // 6 x 7
		                                                    { {0x877788887888ull, 0xd2a5d555d2a5ull}, {0x87777888e088ull, 0xad5a52a58855ull}, {0x788877778776ull, 0xd2a5d555d2a5ull}, {0x788887771f76ull, 0xad5a52a58855ull} } }, // 6 x 8
                                                    { { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 7 x 0
		                                                    { {0x0ull, 0x7full}, {0x0ull, 0x7full}, {0x7eull, 0x7full}, {0x7eull, 0x7full} }, // 7 x 1
		                                                    { {0x3f00ull, 0x3fbfull}, {0x0ull, 0x1ull}, {0xfeull, 0x3fbfull}, {0x0ull, 0x1ull} }, // 7 x 2
		                                                    { {0x1be208ull, 0x1d72baull}, {0x188c30ull, 0x135924ull}, {0x41df6ull, 0x1d72baull}, {0x773ceull, 0x135924ull} }, // 7 x 3
		                                                    { {0x77c3088ull, 0xaaedbaaull}, {0x0ull, 0x1ull}, {0x883cf76ull, 0xaaedbaaull}, {0x0ull, 0x1ull} }, // 7 x 4
		                                                    { {0x1e7e23088ull, 0x492729baaull}, {0x1f88c30e0ull, 0x763592477ull}, {0x6181dcf76ull, 0x492729baaull}, {0x60773cf1eull, 0x763592477ull} }, // 7 x 5
		                                                    { {0x21ddf0c3088ull, 0x34a97b6dbaaull}, {0x0ull, 0x1ull}, {0x1e220f3cf76ull, 0x34a97b6dbaaull}, {0x0ull, 0x1ull} }, // 7 x 6
		                                                    { {0x79f88c307888ull, 0x1249ca6dbd2a5ull}, {0x79f88c307888ull, 0x1db6359242d5aull}, {0x1860773cf8776ull, 0x1249ca6dbd2a5ull}, {0x1860773cf8776ull, 0x1db6359242d5aull} }, // 7 x 7
		                                                    { {0x8777c30c307888ull, 0xd2a5edb6dbd2a5ull}, {0x0ull, 0x1ull}, {0x78883cf3cf8776ull, 0xd2a5edb6dbd2a5ull}, {0x0ull, 0x1ull} } }, // 7 x 8
                                                    { { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 8 x 0
		                                                    { {0x0ull, 0xffull}, {0x0ull, 0xffull}, {0xfeull, 0xffull}, {0xfeull, 0xffull} }, // 8 x 1
		                                                    { {0xfe00ull, 0xff7full}, {0x8888ull, 0xaaaaull}, {0x1feull, 0xff7full}, {0x7776ull, 0xaaaaull} }, // 8 x 2
		                                                    { {0x0ull, 0x1ull}, {0xc30c30ull, 0x924924ull}, {0x0ull, 0x1ull}, {0x3cf3ceull, 0x924924ull} }, // 8 x 3
		                                                    { {0x77788888ull, 0xaaad5aaaull}, {0x78887888ull, 0xd2a5d2a5ull}, {0x88877776ull, 0xaaad5aaaull}, {0x87778776ull, 0xd2a5d2a5ull} }, // 8 x 4
		                                                    { {0x0ull, 0x1ull}, {0x7c3087c308ull, 0xd1245d1245ull}, {0x0ull, 0x1ull}, {0x83cf783cf6ull, 0xd1245d1245ull} }, // 8 x 5
		                                                    { {0x87777888e088ull, 0xd2a5ad5a77aaull}, {0x877788887888ull, 0xad5a2aaa2d5aull}, {0x788887771f76ull, 0xd2a5ad5a77aaull}, {0x788877778776ull, 0xad5a2aaa2d5aull} }, // 8 x 6
		                                                    { {0x0ull, 0x1ull}, {0x8777c30c307888ull, 0xad5a1249242d5aull}, {0x0ull, 0x1ull}, {0x78883cf3cf8776ull, 0xad5a1249242d5aull} }, // 8 x 7
		                                                    { {0x8777788878887888ull, 0xd2a5ad5a2d5ad2a5ull}, {0x8777788878887888ull, 0xad5a52a5d2a52d5aull}, {0x7888877787778776ull, 0xd2a5ad5a2d5ad2a5ull}, {0x7888877787778776ull, 0xad5a52a5d2a52d5aull} } }// 8 x 8
} ;