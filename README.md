================================================================================
【 ソフト名 】ChebyshevPolynomialBdG
【  作成者  】@dc1394
================================================================================

★これは何？
　Bogoliubov-de Gennes方程式を数値的に解くコードです。
　これは、cometscome様の「ChebyshevPolynomialBdG」（
　https://github.com/cometscome/ChebyshevPolynomialBdG ）というコードを、単純に
　C++に移植したものです。以下は、cometscome様による「ChebyshevPolynomialBdG」の
　説明です。

--------------------------------------------------------------------------------
　This solves the Bogoliubov-de Gennes equations and gap equations in the
　s-wave superconductor with the use of the Chebyshev polynomial method. See, 
　Y. Nagai, Y. Ota and M. Machida [arXiv:1105.4939 or
　DOI:10.1143/JPSJ.81.024710] 
　http://journals.jps.jp/doi/abs/10.1143/JPSJ.81.024710
--------------------------------------------------------------------------------

　ビルドには、以下のライブラリが必要です。
　・Eigen

★更新履歴
　2016/08/25 ver.0.1   公開。

★ライセンス
　このソフトはフリーソフトウェアです（2条項BSDライセンス）。
--------------------------------------------------------------------------------
　ChebyshevPolynomialBdG
　Copyright (C) 2016 @dc1394

　ソースコード形式であれバイナリ形式であれ、変更の有無に関わらず、以下の条件を満
　たす限りにおいて、再配布および利用を許可します。

　1.ソースコード形式で再配布する場合、上記著作権表示、本条件書および第2項の責任
　限定規定を必ず含めてください。
　2.バイナリ形式で再配布する場合、上記著作権表示、 本条件書および下記責任限定規
　定を、配布物とともに提供される文書 および/または 他の資料に必ず含めてください。

　本ソフトウェアは著作権者およびコントリビューターによって「現状のまま」提供され
　ており、明示黙示を問わず、商用品として通常そなえるべき品質をそなえているとの保
　証も、特定の目的に適合するとの保証を含め、何の保証もなされません。著作権者もコ
　ントリビューターも、事由のいかんを問わず、損害発生の原因いかんを問わず、かつ責
　任の根拠が契約であるか厳格責任であるか (過失その他)不法行為であるかを問わず、
　仮にそのような損害が発生する可能性を知らされていたとしても、本ソフトウェアの使
　用によって発生した (代替品または代替サービスの提供、使用機会の喪失、データの喪
　失、利益の損失、業務の中断、またそれに限定されない)直接損害、間接損害、偶発的
　な損害、特別損害、懲罰的損害または結果損害のいずれに対しても一切責任を負いませ
　ん。

  Copyright (c) 2016, @dc1394
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
  
  1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer. 
  2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
--------------------------------------------------------------------------------

　ChebyshevPolynomialBdGにはEigen projectによるEigenを使用しています。こちらのラ　イセンスはMPL2になります。
